// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#include "Hardware.Focuser.DIYFocuser.h"

#include <QtSerialPort/QSerialPortInfo>
#include <QDebug>

std::shared_ptr<Focuser> DIYFocuser::Open()
{
    QString portName;
    foreach( const QSerialPortInfo &serialPortInfo, QSerialPortInfo::availablePorts() ) {
        QString desc = serialPortInfo.description();
        if( desc == "USB-SERIAL CH340" || desc == "USB2.0-Serial" || desc == "USB Serial" ) {
            portName = serialPortInfo.portName();
            break;
        }
    }
    if( portName.length() > 0 ) {
        auto focuser = std::make_shared<DIYFocuser>();
        if( focuser->connect( portName ) ) {
            return focuser;
        }
    }
    return 0;
}

bool DIYFocuser::connect( const QString& portName )
{
    serial = new QSerialPort( this );
    serial->setPortName( portName );
    serial->setBaudRate( QSerialPort::Baud9600 );
    serial->setDataBits( QSerialPort::Data8 );
    serial->setParity( QSerialPort::NoParity );
    serial->setStopBits( QSerialPort::OneStop );
    serial->setFlowControl( QSerialPort::NoFlowControl) ;
    if( serial->open( QIODevice::ReadWrite ) ) {
        assert( serial->isOpen() );

        QObject::connect( serial, &QSerialPort::readyRead, this, &DIYFocuser::readSerial ) ;
        qDebug() << "Connected to " << portName;
        return true;
    }
    qDebug() << "Failed to connect to " << portName;
    return false;
}

void DIYFocuser::Close()
{
    if( serial != nullptr && serial->isOpen() ) {
        cancelMoveTo();
        serial->close();
        delete serial;
    }
}

void DIYFocuser::Forward()
{
    if( isInsideMoveTo() ) {
        cancelMoveTo();
        return;
    }
    focuserPos = INT_MIN;
    writeToSerial( QString( "FWD %1\n" ).arg( QString::number( stepsToGo ) ) );
}

void DIYFocuser::Backward()
{
    if( isInsideMoveTo() ) {
        cancelMoveTo();
        return;
    }
    focuserPos = INT_MIN;
    writeToSerial( QString( "BWD %1\n" ).arg( QString::number( stepsToGo ) ) );
}

void DIYFocuser::MarkZero()
{
    cancelMoveTo();
    writeToSerial( QString( "POS ZERO\n" ) );
}

void DIYFocuser::GoToPos( int pos )
{
    if( !isInsideMoveTo() ) {
        focuserPos = INT_MIN;
        writeToSerial( QString( "POS GET\n" ) );
    }
    targetPos = pos;
}

void DIYFocuser::writeToSerial( const QString& str ) const
{
    assert( serial->isWritable() );
    serial->write( str.toLocal8Bit().constData() );
    serial->waitForBytesWritten( 5000) ;
}

void DIYFocuser::readSerial() const
{
    QByteArray data = serial->readAll();
    QString reply = QString::fromStdString( data.toStdString() );
    qDebug() << reply;
    if( reply.contains( "RESET" ) ) {
        writeToSerial( QString( "POS GET\n" ) );
        return;
    }
    if( reply.startsWith( "POS " ) ) {
        int pos = reply.indexOf( ' ' );
        if( pos != -1 ) {
            focuserPos = reply.mid( pos + 1 ).toInt();
        }
    }
    if( isInsideMoveTo() ) {
        if( focuserPos < targetPos ) {
            writeToSerial( QString( "FWD %1\n" ).arg( QString::number( qMin( 128, targetPos - focuserPos ) ) ) );
        } else if( focuserPos > targetPos ) {
            writeToSerial( QString( "BWD %1\n" ).arg( QString::number( qMin( 128, focuserPos - targetPos ) ) ) );
        } else {
            cancelMoveTo();
        }
    }
}
