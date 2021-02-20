// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#ifndef IMAGE_QT_H
#define IMAGE_QT_H

#include <QImage>
#include <QPixmap>

#include "Image.Formats.h"
#include "Image.Math.h"

#include <Image.RgbImage.h>

inline QImage qImage( const CRgbImage* image )
{
    return QImage( image->RgbPixels(), image->Width(), image->Height(), image->ByteWidth(), QImage::Format_RGB888 );
}

inline QPixmap qPixmap( const CRgbImage* image )
{
    return QPixmap::fromImage( QImage( image->RgbPixels(), image->Width(), image->Height(), image->ByteWidth(), QImage::Format_RGB888 ) );
}

inline QPixmap qPixmap( std::shared_ptr<const CRgbImage> image )
{
    return qPixmap( image.get() );
}

#endif // IMAGE_QT_H