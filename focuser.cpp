// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#include "focuser.h"

#include <QtSerialPort/QSerialPortInfo>
#include <QDebug>

bool Focuser::Open()
{
    QString portName;
    foreach( const QSerialPortInfo &serialPortInfo, QSerialPortInfo::availablePorts() ) {
        QString desc = serialPortInfo.description();
        if( desc == "USB-SERIAL CH340" || desc == "USB2.0-Serial" ) {
            portName = serialPortInfo.portName();
            break;
        }
    }
    if( portName.length() > 0 ) {
        serial = new QSerialPort( this );
        serial->setPortName( portName );
        serial->setBaudRate( QSerialPort::Baud9600 );
        serial->setDataBits( QSerialPort::Data8 );
        serial->setParity( QSerialPort::NoParity );
        serial->setStopBits( QSerialPort::OneStop );
        serial->setFlowControl( QSerialPort::NoFlowControl) ;
        serial->open( QIODevice::ReadWrite );
        assert( serial->isOpen() );

        QObject::connect( serial, &QSerialPort::readyRead, this, &Focuser::readSerial ) ;

        qDebug() << "Connected to " << portName;
        return true;
    }
    return false;
}

void Focuser::Close()
{
    if( serial != nullptr && serial->isOpen() ) {
        serial->close();
        delete serial;
    }
}

void Focuser::Forward()
{
    writeToSerial( QString( "FWD %1\n" ).arg( QString::number( stepsToGo ) ) );
}

void Focuser::Backward()
{
    writeToSerial( QString( "BWD %1\n" ).arg( QString::number( stepsToGo ) ) );
}

void Focuser::writeToSerial( const QString& str )
{
    assert( serial->isWritable() );
    serial->write( str.toLocal8Bit().constData() );
    serial->waitForBytesWritten( 5000) ;
}

void Focuser::readSerial()
{
    QByteArray data = serial->readAll();
    qDebug() << QString::fromStdString( data.toStdString() );
}
