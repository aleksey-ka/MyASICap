// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#ifndef MAINFRAME_H
#define MAINFRAME_H

#include <QMainWindow>

#include <QFutureWatcher>
#include <QTimer>
#include <QSettings>

#include "camera.h"
#include "focuser.h"

#include "focusercontrols.h"

namespace Ui {
    class MainFrame;
}

class MainFrame : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainFrame( QWidget *parent = nullptr );
    ~MainFrame();

private slots:
    void on_closeButton_clicked();

    void on_toggleFullScreenButton_clicked();

    void on_cameraSelectionCombo_currentIndexChanged( int );

    void on_cameraOpenCloseButton_clicked();

    void on_exposureSpinBox_valueChanged( int value );

    void on_exposureSlider_valueChanged( int value );

    void on_gainSpinBox_valueChanged( int value );

    void on_gainSlider_valueChanged( int value );

    void on_captureButton_clicked();

    void on_guideUp();

    void on_guideDown();

    void on_guideLeft();

    void on_guideRight();

    void on_coolerCheckBox_stateChanged( int );

    void on_temperatureSpinBox_valueChanged( int );

private:
    Ui::MainFrame *ui;
    void updateUI();

    QSettings settings;

    std::vector<std::shared_ptr<ASI_CAMERA_INFO>> camerasInfo;
    std::shared_ptr<Camera> camera;

    std::shared_ptr<ASI_CAMERA_INFO> openCamera( int index );
    void closeCamera();

    // Capture
    QFutureWatcher<std::shared_ptr<const Raw16Image>> imageReadyWatcher;
    QFutureWatcher<QString> imageSavedWatcher;
    int exposureRemainingTime;
    int capturedFrames = 0;
    QTimer exposureTimer;
    void startCapture();
    void showCaptureStatus();
    void imageReady();
    QString saveToPath;
    void imageSaved();

    // Manual guider controls
    int guiding = -1;
    void guide( ASI_GUIDE_DIRECTION );
    void guideStop();

    // Focuser
    Focuser focuser;

    // Rendering
    ulong render( const ushort* raw, int width, int height );
    bool drawTargetingCircle = false;

    // Exposure controls and scaling
    void setExposureInSpinBox( int exposure );
    static int sliderPosToExposure( int pos );
    static int exposureToSliderPos( int exposure );
    static int exposureSuffixToScale( const QString& suffix );
    static int exposureToScaleAndSuffix( int exposure, QString& suffix );
    static QString exposureToString( int exposure );
    static int exposureFromValueAndSuffix( int value, const QString& suffix );

    // Temperature control
    void setCooler( bool isOn, int targetTemperature );
};

#endif // MAINFRAME_H
