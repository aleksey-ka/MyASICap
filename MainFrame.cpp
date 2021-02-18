// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#include "MainFrame.h"
#include "ui_MainFrame.h"
#include <QShortcut>
#include <QtConcurrent/QtConcurrent>
#include <QPainter>
#include <QInputDialog>
#include <QDebug>

#include "Camera.ZWO.h"
#include "Camera.Mock.h"
#include "Image.Formats.h"
#include "Renderer.h"

#include <chrono>

MainFrame::MainFrame( QWidget *parent ) :
    QMainWindow( parent ),
    ui( new Ui::MainFrame )
{
    // Multithreading noticibly improves throughput especially when writing to compressed image formats
    int numberOfCores = QThread::idealThreadCount();
    QThreadPool::globalInstance()->setMaxThreadCount( std::min( numberOfCores, 4 ) );

    setWindowIcon( QIcon( ":MainFrame.ico" ) );

    ui->setupUi( this ) ;

    ui->closeButton->hide();

    ui->imageView->setBackgroundRole( QPalette::Base );
    ui->imageView->setSizePolicy( QSizePolicy::Ignored, QSizePolicy::Ignored );
    ui->imageView->setScaledContents( false );
    ui->imageView->setAlignment( Qt::AlignCenter );

    ui->histogramView->setBackgroundRole( QPalette::Base );
    ui->histogramView->setSizePolicy( QSizePolicy::Ignored, QSizePolicy::Ignored );
    ui->histogramView->setScaledContents( false );
    ui->histogramView->setAlignment( Qt::AlignCenter );

    new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_F ), this, SLOT( on_toggleFullScreenButton_clicked() ) );
    new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_Q ), this, SLOT( on_closeButton_clicked() ) );

    new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_Up ), this, SLOT( on_guideUp() ) );
    new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_Down ), this, SLOT( on_guideDown() ) );
    new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_Left ), this, SLOT( on_guideLeft() ) );
    new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_Right ), this, SLOT( on_guideRight() ) );

    // TODO: In Qt 5.15 lambdas can be used in QShortcut constructor
    connect( new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_T ), this ), &QShortcut::activated, [=]() { drawTargetingCircle = !drawTargetingCircle; } );

    connect( &imageReadyWatcher, &QFutureWatcher<std::shared_ptr<Raw16Image>>::finished, this, &MainFrame::imageReady );
    connect( &imageSavedWatcher, &QFutureWatcher<QString>::finished, this, &MainFrame::imageSaved );

    connect( &exposureTimer, &QTimer::timeout, [=]() { if( exposureRemainingTime > 0 ) exposureRemainingTime--; showCaptureStatus(); } );

    int count = ASICamera::GetCount();
    for( int i = 0; i < count; i++ ) {
        camerasInfo.emplace_back( ASICamera::GetInfo( i ) );
    }

    count = MockCamera::GetCount();
    for( int i = 0; i < count; i++ ) {
        camerasInfo.emplace_back( MockCamera::GetInfo( i ) );
    }

    // TODO: Fixing a bug with text color on Raspberry Pi (old Qt?). It shows always gray
    // To fix it needs changing the combo to editable and the edit inside the combo to read-only
    ui->cameraSelectionCombo->lineEdit()->setReadOnly( true );

    for( size_t i = 0; i < camerasInfo.size(); i++ ) {
        ui->cameraSelectionCombo->addItem( camerasInfo[i]->Name, QVariant( i ) );
    }

    ui->objectNameEdit->setText( settings.value( "Name" ).toString() );
    ui->formatComboBox->setCurrentText( settings.value( "FileFormat", ".png" ).toString() );

    auto defaultPath = QStandardPaths::writableLocation( QStandardPaths::DocumentsLocation ) + QDir::separator() +
            QApplication::applicationName() + QDir::separator() + "{TIME}{NAME}{FILTER}";
    saveToPath = settings.value( "SaveTo", defaultPath ).toString();
    ui->saveToEdit->setText( saveToPath );
    ui->saveToEdit->home( false );

    setExposureInSpinBox( settings.value( "Exposure", 100000 ).toInt() );
    ui->gainSpinBox->setValue( settings.value( "Gain", 0 ).toInt() );
    ui->offsetSpinBox->setValue( settings.value( "Offset", 64 ).toInt() );
    ui->useCameraWhiteBalanceCheckBox->setChecked( settings.value( "UseCameraWhiteBalance", false ).toBool() );

   if( focuser.Open() ) {
       // TODO: In Qt 5.15 lambdas can be used in QShortcut constructor
       connect( new QShortcut( QKeySequence( Qt::CTRL + Qt::SHIFT + Qt::Key_Up ), this ), &QShortcut::activated, [=]() { focuser.Forward(); } );
       connect( new QShortcut( QKeySequence( Qt::CTRL + Qt::SHIFT + Qt::Key_Down ), this ), &QShortcut::activated, [=]() { focuser.Backward(); } );
       connect( new QShortcut( QKeySequence( Qt::CTRL + Qt::SHIFT + Qt::Key_Right ), this ), &QShortcut::activated, [=]() { focuser.StepUp(); } );
       connect( new QShortcut( QKeySequence( Qt::CTRL + Qt::SHIFT + Qt::Key_Left ), this ), &QShortcut::activated, [=]() { focuser.StepDown(); } );
       connect( new QShortcut( QKeySequence( Qt::CTRL + Qt::SHIFT + Qt::Key_Z ), this ), &QShortcut::activated, [=]() { focuser.MarkZero(); } );
       connect( new QShortcut( QKeySequence( Qt::CTRL + Qt::SHIFT + Qt::Key_0 ), this ), &QShortcut::activated, [=]() { focuser.GoToPos( 0 ); } );
       connect( new QShortcut( QKeySequence( Qt::CTRL + Qt::SHIFT + Qt::Key_G ), this ),
            &QShortcut::activated, [=]() {
                bool ok = false;
                int pos = QInputDialog::getInt( this, "Focuser", "MoveTo:", 0, INT_MIN, INT_MAX, 1, &ok );
                if( ok ) {
                    focuser.GoToPos( pos );
                }
            }
       );
   }

   if( filterWheel.Open() ) {
       // TODO: Fixing a bug with text color on Raspberry Pi (old Qt?). It shows always gray
       // To fix it needs changing the combo to editable and the edit inside the combo to read-only
       ui->filterWheelComboBox->lineEdit()->setReadOnly( true );
       ui->filterComboBox->lineEdit()->setReadOnly( true );

       auto filterWheelDefs = QDir( "", "*.FilterWheel" ).entryList( QDir::Files );
       if( filterWheelDefs.size() > 0 ) {
           QSignalBlocker lock( ui->filterWheelComboBox );
           int index = 0;
            for( size_t i = 0; i < filterWheelDefs.size(); i++ ) {
                auto name = QFileInfo( filterWheelDefs[i] ).baseName();
                ui->filterWheelComboBox->addItem( name );
                if( name == settings.value( "FilterWheel" ) ) {
                    index = i;
                }
            }
            ui->filterWheelComboBox->setCurrentIndex( index );
            on_filterWheelComboBox_currentIndexChanged( index );
       } else {
            ui->filterWheelComboBox->addItem( "Filter Wheel" );
       }

       connect( ui->filterComboBox, QOverload<int>::of( &QComboBox::currentIndexChanged ),
            [this]( int index ) { filterWheel.SetPosition( index ); } );
   }

   ui->imageView->setCursor( QCursor( QPixmap( ":CrossHair.png" ), 23, 23 ) );

   updateUI();
}

MainFrame::~MainFrame()
{
    if( camera != 0 ) {
        closeCamera();
    }
    focuser.Close();
    filterWheel.Close();
    delete ui;
}

void MainFrame::updateUI()
{
    ui->histogramView->setVisible( camera != 0 );
    ui->cameraControlsFrame->setVisible( camera != 0 );
    ui->captureFrame->setVisible( camera != 0 );
    ui->temperatureFrame->setVisible( camera != 0 && camera->HasCooler() );
    ui->cameraOpenCloseButton->setText( camera != 0 ? "X" : ">" );
    ui->useCameraWhiteBalanceCheckBox->setVisible( camera != 0 && camera->GetInfo()->IsColorCam );
    ui->filterWheelFrame->setVisible( filterWheel.GetSlotsCount() > 0 );
}

void MainFrame::on_closeButton_clicked()
{
    close();
}

void MainFrame::on_toggleFullScreenButton_clicked()
{
    bool isFullScreen = windowState().testFlag( Qt::WindowFullScreen );
    if( isFullScreen ) {
        ui->closeButton->hide();
        showNormal();  
    } else {
        showFullScreen();
        ui->closeButton->show();
    }
}

void MainFrame::on_cameraSelectionCombo_currentIndexChanged( int index )
{
    if( camera != 0 ) {
        closeCamera();
        updateUI();
    }

    auto info = camerasInfo[index];

    const static QString namedValue( "%1: <span style='color:#008800;'>%2</span><br>");

    QString txt;
    txt.append( namedValue.arg( "Bit depth", QString::number( info->BitDepth ) ) );
    txt.append( namedValue.arg( "Pixel size", QString::number( info->PixelSize ) ) );

    ui->infoLabel->setText( txt );

    ui->imageView->clear();
    ui->imageView->setText( "No image" );
}

void MainFrame::on_cameraOpenCloseButton_clicked()
{
    if( camera == 0 ) {
        int selectedIndex = ui->cameraSelectionCombo->currentData().toInt();
        auto cameraInfo = openCamera( selectedIndex );
        setWindowTitle( cameraInfo->Name );
    } else {
        closeCamera();
        setWindowTitle( "No camera" );
    }
    updateUI();
}

void MainFrame::on_exposureSpinBox_valueChanged( int value )
{
    int exposure = exposureFromValueAndSuffix( value, ui->exposureSpinBox->suffix() );
    ui->exposureSlider->setValue( exposureToSliderPos( exposure ) );
}

void MainFrame::on_exposureSlider_valueChanged( int pos )
{
    int exposure = sliderPosToExposure( pos );
    setExposureInSpinBox( exposure );
}

void MainFrame::on_gainSpinBox_valueChanged( int value )
{
    ui->gainSlider->setValue( value );
}

void MainFrame::on_gainSlider_valueChanged( int value )
{
    ui->gainSpinBox->setValue( value );
}

void MainFrame::on_coolerCheckBox_stateChanged( int state )
{
    setCooler( state == Qt::Checked, ui->temperatureSpinBox->value() );
}

void MainFrame::on_temperatureSpinBox_valueChanged( int targetTemperature )
{
    setCooler( ui->coolerCheckBox->isChecked(), targetTemperature );
}

void MainFrame::on_guideUp()
{
    guide( ASI_GUIDE_NORTH );
}

void MainFrame::on_guideDown()
{
    guide( ASI_GUIDE_SOUTH );
}

void MainFrame::on_guideLeft()
{
    guide( ASI_GUIDE_EAST );
}

void MainFrame::on_guideRight()
{
    guide( ASI_GUIDE_WEST );
}

void MainFrame::guide( ASI_GUIDE_DIRECTION direction )
{
    if( guiding != -1 ) {
        camera->GuideOff( (ASI_GUIDE_DIRECTION)guiding );
    }
    if( guiding != direction ) {
        camera->GuideOn( direction );
        guiding = direction;
    } else {
        guiding = -1;
    }
}

void MainFrame::guideStop()
{
    if( guiding != -1 ) {
        camera->GuideOff( (ASI_GUIDE_DIRECTION)guiding );
        guiding = -1;
    }
}

void MainFrame::on_captureButton_clicked()
{
    ui->captureButton->setEnabled( false );
    ui->saveToCheckBox->setEnabled( false );

    if( ui->saveToCheckBox->isChecked() ) {
        auto name = ui->objectNameEdit->text();
        settings.setValue( "Name", name );
        auto path = ui->saveToEdit->text();
        settings.setValue( "SaveTo", path );

        static const QString TIMESTAMP( "{TIME}" );
        auto timestamp = QDateTime::currentDateTime().toString("yyyy-MM-ddThh-mm-ss");
        if( path.contains( TIMESTAMP ) ) {
            path.replace( TIMESTAMP, "-" + timestamp );
        }
        static const QString FILTER( "{FILTER}" );
        if( path.contains( FILTER ) ) {
            auto filter = ui->filterComboBox->currentText();
            if( not filter.isEmpty() ) {
                path.replace( FILTER, "-" + filter );
            } else {
                path.replace( FILTER, "" );
            }
        }
        static const QString NAME( "{NAME}" );
        if( path.contains( NAME ) ) {
            if( not name.isEmpty() ) {
                path.replace( NAME, "-" + name );
            } else {
                path.replace( NAME, "" );
            }
        }
        path.replace( "\\-", "\\" );
        path.replace( "/-", "/" );
        saveToPath = path;
    } else {
        saveToPath.clear();
    }
    seriesId = std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() );

    settings.setValue( "Exposure", ui->exposureSpinBox->value() * exposureSuffixToScale( ui->exposureSpinBox->suffix() ) );
    settings.setValue( "Gain", ui->gainSpinBox->value() );
    settings.setValue( "Offset", ui->offsetSpinBox->value() );
    settings.setValue( "UseCameraWhiteBalance", ui->useCameraWhiteBalanceCheckBox->isChecked() );

    startCapture();
}

std::shared_ptr<ASI_CAMERA_INFO> MainFrame::openCamera( int index )
{
    auto start = std::chrono::steady_clock::now();

    auto cameraInfo = camerasInfo[index];
    if( cameraInfo->CameraID < 0 ) {
        // Mock camera
        camera = MockCamera::Open( cameraInfo->CameraID );
    } else {
        camera = ASICamera::Open( cameraInfo->CameraID );
    }

    int width = 0;
    int height = 0;
    int bin = 0;
    ASI_IMG_TYPE imgType = ASI_IMG_END;
    camera->GetROIFormat( width, height, bin, imgType );
    imgType = ASI_IMG_RAW16;
    camera->SetROIFormat( width, height, bin, imgType );
    camera->GetROIFormat( width, height, bin, imgType );
    qDebug() << width << "x" << height << " bin" << bin;
    switch( imgType ) {
        case ASI_IMG_RAW8: qDebug() << "RAW8"; break;
        case ASI_IMG_RGB24: qDebug() << "RGB24"; break;
        case ASI_IMG_RAW16: qDebug() << "RAW16"; break;
        case ASI_IMG_Y8: qDebug() << "Y8"; break;
        default:
            assert( false );
    }

    //camera->PrintDebugInfo();

    auto end = std::chrono::steady_clock::now();
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>( end - start ).count();
    qDebug() << "Camera initialized in" << msec << "msec";

    return cameraInfo;
}

void MainFrame::closeCamera()
{
    guideStop();
    camera->Close();
    camera = 0;
}

void MainFrame::startCapture()
{
    auto exposure = exposureFromValueAndSuffix( ui->exposureSpinBox->value(), ui->exposureSpinBox->suffix() );
    auto gain = ui->gainSpinBox->value();
    auto offset = ui->offsetSpinBox->value();
    auto useCameraWhiteBalance = ui->useCameraWhiteBalanceCheckBox->isChecked();

    camera->SetExposure( exposure );
    camera->SetGain( gain );
    camera->SetWhiteBalanceR( useCameraWhiteBalance ? 52 : 50 );
    camera->SetWhiteBalanceB( useCameraWhiteBalance ? 95 : 50 );
    camera->SetOffset( offset );

    if( filterWheel.GetSlotsCount() > 0 ) {
        auto channel = ui->filterComboBox->currentText();
        camera->SetChannel( channel.toStdString().c_str() );
        auto fullFilterDescription = ui->filterComboBox->currentData( Qt::ToolTipRole ).toString();
        camera->SetFilterDescription( fullFilterDescription.toStdString().c_str() );
    }
    camera->SetSeriesId( seriesId );

    auto info = camera->GetInfo();
    ui->ePerADULabel->setText( QString( "e<sup>-</sup>/ADU: <span style='color:#008800;'>%1</span>" ).arg( QString::number( info->ElecPerADU, 'f', 3 ) ) );

    imageReadyWatcher.setFuture( QtConcurrent::run( [=]() {
        auto start = std::chrono::steady_clock::now();

        auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>( start.time_since_epoch() ).count();
        if( capturedFrames > 0 ) {
            qDebug() << "Since previous exposure " << timestamp - previousTimestamp << "msec";
            qDebug() << "FPS " << ( 1000.0 * capturedFrames ) / ( timestamp - startTimestamp );
        } else {
            startTimestamp = timestamp;
        }
        previousTimestamp = timestamp;

        auto result = camera->DoExposure();

        auto msec = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() - start ).count();
        qDebug() << "Captured in " << msec << "msec";

        return result;
    } ) );

    exposureRemainingTime = exposure / 1000000;
    exposureTimer.start( 1000 );
    showCaptureStatus();
}

void MainFrame::showCaptureStatus()
{
    ui->captureButton->setText( QString::number( exposureRemainingTime ) +
        ( capturedFrames > 0 ? " (" + QString::number( capturedFrames ) + ")" : "" ) );
}

void MainFrame::imageReady()
{
    exposureTimer.stop();

    auto result = imageReadyWatcher.result();

    int currentIndex = capturedFrames.fetchAndAddOrdered( 1 );

    if( result != 0 ) {
        imageSavedWatcher.setFuture( QtConcurrent::run( [=]() {
            const static QString namedValue( "%1: <span style='color:#008800;'>%2</span>%3<br>");

            QString txt;
            const auto& info = result->Info();
            txt.append( namedValue.arg( "Size", QString::number( info.Width ) + "x" + QString::number( info.Height ), "" ) );
            txt.append( namedValue.arg( "Exposure", exposureToString( info.Exposure ), "" ) );
            txt.append( namedValue.arg( "Gain", QString::number( info.Gain ), "" ) );
            txt.append( namedValue.arg( "Offset", QString::number( info.Offset ), "" ) );
            txt.append( namedValue.arg( "Temperature", QString::number( info.Temperature, 'f', 1 ), " <sup>0</sup>C" ) );
            if( not info.FilterDescription.empty() ) {
                txt.append( namedValue.arg( "Filter", QString( info.FilterDescription.c_str() ), "" ) );
            }
            txt.append( "<br>" );

            if( saveToPath.length() > 0 ) {
                auto ext = ui->formatComboBox->currentText();
                settings.setValue( "FileFormat", ext );
                const ImageFileFormat* format = 0;
                if( ext == ".png") {
                    const static Png16BitGrayscale png16BitGrayscale;
                    format = &png16BitGrayscale;
                } else {
                    assert( ext == ".pixels" );
                }

                QDir().mkpath( saveToPath );
                const static QString nameTemplate( "%1.%2%3%4.u16%5" );
                auto name = nameTemplate.arg( QString::number( info.SeriesId, 16 ), QString::number( currentIndex ).rightJustified( 5, '0' ),
                    ( info.CFA.empty() ? "" : ".cfa" ), ( info.Channel.empty() ? "" : "." + info.Channel ).c_str(), ext );

                result->SaveToFile( ( saveToPath + QDir::separator() + name ).toStdString().c_str(), format );
                txt.append( namedValue.arg( "Saved As", name, "" ) );
            }

            return txt;

        } ) );

        auto msec = render( result->RawPixels(), result->Width(), result->Height() );
        qDebug() << "Rendered in " << msec << "msec";
    }

    if( ui->continuousCaptureCheckBox->isChecked() ) {
        ui->continuousCaptureCheckBox->setText( "Continuous capture (uncheck to stop)" );
        startCapture();
    } else {
        capturedFrames = 0;
        ui->continuousCaptureCheckBox->setText( "Continuous capture" );
        ui->captureButton->setText( "Capture" );
        ui->captureButton->setEnabled( true );
        ui->saveToCheckBox->setEnabled( true );
    }
}

void MainFrame::imageSaved()
{
    ui->infoLabel->setText( imageSavedWatcher.result() );
}

ulong MainFrame::render( const ushort* raw, int width, int height )
{
    auto start = std::chrono::steady_clock::now();

    Renderer renderer( raw, width, height );
    QPixmap pixmap = renderer.Render( ui->showFullResolution->isChecked() ? RM_FullResolution : RM_HalfResolution );
    if( drawTargetingCircle ) {
       QPainter painter( &pixmap );
       QPen pen( QColor::fromRgb( 0xFF, 0, 0 ) );
       pen.setWidth( 3 );
       painter.setPen( pen );
       painter.drawEllipse( pixmap.width() / 2 - 30, pixmap.height() / 2 - 30, 60, 60 );
    }
    ui->imageView->setPixmap( pixmap );
    ui->histogramView->setPixmap( renderer.RenderHistogram() );

    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>( end - start ).count();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Exposure scaling

void MainFrame::setExposureInSpinBox( int exposure )
{
    // First set the slider to the appropriate range
    ui->exposureSlider->setValue( exposureToSliderPos( exposure ) );

    // Now set the exposure in the exposure spin box to the exact value
    // with convinient time units (s/ms/us)
    QString suffix;
    exposure /= exposureToScaleAndSuffix( exposure, suffix );
    ui->exposureSpinBox->setSuffix( suffix );
    ui->exposureSpinBox->setValue( exposure );
}

// Slider position to exposure translation table (quasi logorithmic)
static int exposures[] =
    { 1, 1, 1, 2, 2, 3, 4, 5, 6, 8,
      10, 12, 15, 20, 25, 30, 40, 50, 60, 80,
      100, 120, 150, 200, 250, 300, 400, 500, 600, 800,
      1000, 1200, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000,
      10000, 12000, 15000, 20000, 25000, 30000, 40000, 50000, 60000, 80000,
      100000, 120000, 150000, 200000, 250000, 300000, 400000, 500000, 600000, 800000,
      1000000, 1200000, 1500000, 2000000, 2500000, 3000000, 4000000, 5000000, 6000000, 8000000,
      10000000, 12000000, 15000000, 20000000, 25000000, 30000000, 40000000, 50000000, 60000000, 80000000,
      100000000 };

int MainFrame::sliderPosToExposure( int pos )
{
    return exposures[pos];
}

int MainFrame::exposureToSliderPos( int exposure )
{
    const int lastPos = sizeof( exposures ) / sizeof( int ) - 1;
    for( int i = 0; i <= lastPos; i++ ) {
        if( exposures[i] > exposure ) {
            return i - 1;
        }
    }
    return lastPos;
}

int MainFrame::exposureSuffixToScale( const QString& suffix )
{
    if( suffix == " s" ) {
        return 1000000;
    } else if( suffix == " ms" ) {
        return 1000;
    } else {
        assert( suffix == " us");
        return 1;
    }
}

int MainFrame::exposureToScaleAndSuffix( int exposure, QString& suffix )
{
    if( exposure >= 10000000 ) {
        // Seconds for exposures equal or above 10 s
        suffix = " s";
        return 1000000;
    } else if( exposure >= 10000 ) {
        // Milliseconds for exposures equal or above 10 ms
        suffix = " ms";
        return 1000;
    } else {
        // Microseconds for exposures below 10 ms
        suffix = " us";
        return 1;
    }
}

QString MainFrame::exposureToString( int exposure )
{
    QString suffix;
    int scale = exposureToScaleAndSuffix( exposure, suffix );
    return QString::number( exposure / scale ) + suffix;
}

int MainFrame::exposureFromValueAndSuffix( int value, const QString& suffix )
{
    return value *  exposureSuffixToScale( suffix );
}

 // Temperature control
void MainFrame::setCooler( bool coolerOn, int targetTemperature )
{
    if( coolerOn ) {
        camera->SetTargetTemperature( targetTemperature );
    }
    camera->SetCoolerOn( coolerOn );
}

void MainFrame::on_filterWheelComboBox_currentIndexChanged( int index )
{
    auto name = ui->filterWheelComboBox->itemText( index );
    auto fileName = name + ".FilterWheel";

    QStringList slotNames;
    if( QFileInfo::exists( fileName ) ) {
        QFile inputFile( fileName );
        if( inputFile.open( QIODevice::ReadOnly ) ) {
           QTextStream in( &inputFile );
           while( !in.atEnd() ) {
              slotNames.append( in.readLine() );
           }
           inputFile.close();
        }
        settings.setValue( "FilterWheel", name );
    }

    QSignalBlocker lock( ui->filterComboBox );
    ui->filterComboBox->clear();
    for( size_t i = 0; i < filterWheel.GetSlotsCount(); i++ ) {
        QString name;
        QString description;
        if( i < slotNames.size() ) {
            QString slot = slotNames.at( i );
            name = slot.section( " ", 0, 0 );
            description = slot.section( " ", 1 );
        } else {
            name = QString::number( i + 1 );
        }
        ui->filterComboBox->addItem( name );
        ui->filterComboBox->setItemData( i, description, Qt::ToolTipRole );

    }
    ui->filterComboBox->setCurrentIndex( filterWheel.GetPosition() );
}
