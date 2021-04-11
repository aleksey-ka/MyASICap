// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#include "Renderer.h"

#include <Image.Qt.h>
#include <Image.Debayer.HalfRes.h>
#include <Image.Debayer.HQLinear.h>
#include <Image.Debayer.CFA.h>

#include <QPainter>

Renderer::Renderer( const ushort* _raw, int _width, int _height, int _bitDepth ) :
    raw( _raw ), width( _width ), height( _height ), bitDepth( _bitDepth )
{

}

QPixmap Renderer::Render( TRenderingMethod method, int x, int y, int W, int H )
{
    // Initialize histogram
    const int hSize = 256;
    histR.resize( hSize );
    histG.resize( hSize );
    histB.resize( hSize );

    if( method == RM_HalfResolution ) {
        x -= W / 2;
        y -= H / 2;
        size_t w = W > 0 ? W : width / 2;
        size_t h = H > 0 ? H : height / 2;
        size_t byteWidth = 3 * w;
        std::vector<uchar> pixels( byteWidth * h );   
        uchar* rgb = pixels.data();

        CDebayer_RawU16_HalfRes debayer( raw, width, height, bitDepth );
        debayer.ToRgbU8( rgb, byteWidth, x, y, w, h, histR.data(), histG.data(), histB.data() );
        maxValue = debayer.MaxValue;
        maxCount = debayer.MaxCount;
        minValue = debayer.MinValue;
        minCount = debayer.MinCount;

        return Qt::CreatePixmap( rgb, w, h, byteWidth );
    } else if( method == RM_FullResolution ) {
        assert( method == RM_FullResolution );

        size_t w = W > 0 ? W : width;
        size_t h = H > 0 ? H : height;
        size_t byteWidth = 3 * w;
        std::vector<uchar> pixels( byteWidth * h );
        uchar* rgb = pixels.data();

        CDebayer_RawU16_HQLinear debayer( raw, width, height, bitDepth );
        debayer.ToRgbU8( rgb, byteWidth, x, y, w, h, histR.data(), histG.data(), histB.data() );
        maxValue = debayer.MaxValue;
        maxCount = debayer.MaxCount;
        minValue = debayer.MinValue;
        minCount = debayer.MinCount;

        return Qt::CreatePixmap( rgb, w, h, byteWidth );
    } else {
        assert( method == RM_CFA );

        size_t w = W > 0 ? W : width;
        size_t h = H > 0 ? H : height;
        size_t byteWidth = 3 * w;
        std::vector<uchar> pixels( byteWidth * h );
        uchar* rgb = pixels.data();

        CDebayer_RawU16_CFA debayer( raw, width, height, bitDepth );
        debayer.ToRgbU8( rgb, byteWidth, x, y, w, h, histR.data(), histG.data(), histB.data() );
        maxValue = debayer.MaxValue;
        maxCount = debayer.MaxCount;
        minValue = debayer.MinValue;
        minCount = debayer.MinCount;

        return Qt::CreatePixmap( rgb, w, h, byteWidth );
    }
}

QPixmap Renderer::RenderQ( int x, int y, int W, int H )
{
    // Initialize histogram
    const int hSize = 256;
    histR.resize( hSize );
    histG.resize( hSize );
    histB.resize( hSize );

    x -= W / 2;
    y -= H / 2;
    size_t w = W > 0 ? W : width / 2;
    size_t h = H > 0 ? H : height / 2;

    size_t byteWidth = 3 * w;
    std::vector<uchar> pixels( byteWidth * h );
    uchar* rgb = pixels.data();

    CDebayer_RawU16_HalfRes debayer( raw, width, height, bitDepth );
    debayer.ToRgbU8( rgb, byteWidth, x, y, w, h, histR.data(), histG.data(), histB.data() );
    maxValue = debayer.MaxValue;
    maxCount = debayer.MaxCount;
    minValue = debayer.MinValue;
    minCount = debayer.MinCount;

    size_t byteWidth2= byteWidth / 2;
    std::vector<uchar> pixels2( byteWidth2 * h / 2);
    uchar* rgb2 = pixels.data();
    for( int i = 0; i < h / 2; i++ ) {
        const uchar* ptr = rgb + 2 * i * byteWidth;
        uchar* ptr2 = rgb2 + i * byteWidth2;
        for( int j = 0; j < w / 2; j++ ) {
            const uchar* p = ptr + 6 * j;
            uchar* p2 = ptr2 + 3 * j;
            p2[0] = ( p[0] + p[3] + p[byteWidth] + p[byteWidth + 3] ) / 4;
            p2[1] = ( p[1] + p[4] + p[byteWidth + 1] + p[byteWidth + 4] ) / 4;
            p2[2] = ( p[2] + p[5] + p[byteWidth + 2] + p[byteWidth + 5] ) / 4;
        }
    }

    return Qt::CreatePixmap( rgb2, w / 2, h / 2, byteWidth2 );
}


static QString collapseNumber( uint n )
{
    // Collapse uint to Kilo-s and Mega-s with at most 3 digits (1.24K, 24.5M)
    if( n < 1000 ) {
        return QString::number( n );
    } else if( n < 10000 ) {
        return QString::number( n / 1000.0, 'f', 2 ) + "K";
    } else if( n < 100000 ) {
        return QString::number( n / 1000.0, 'f', 1 ) + "K";
    } else if( n < 1000000 ) {
        return QString::number( n / 1000 ) + "K";
    } else if( n < 10000000 ) {
        return QString::number( n / 1000000.0, 'f', 2 ) + "M";
    } else if( n < 100000000 ) {
        return QString::number( n / 1000000.0, 'f', 1 ) + "M";
    } else {
        return QString::number( n / 1000000 ) + "M";
    }
}

QPixmap Renderer::RenderHistogram()
{
    const uint* r = histR.data();
    const uint* g = histG.data();
    const uint* b = histB.data();

    int size = histR.size();

    uint max = 0;
    for( int i = 0; i < size - 1; i++ ) {
        if( r[i] > max ) {
            max = r[i];
        }
        if( g[i] / 2 > max ) {
            max = g[i] / 2;
        }
        if( b[i] > max ) {
            max = b[i];
        }
    }
    if( max == 0 ) {
        max = std::max( std::max( r[size - 1], g[size - 1] ), b[size - 1 ] );
    }

    int h = 128;
    int biteWidth = 3 * size;
    std::vector<uchar> hist( biteWidth * h );
    std::fill( hist.begin(), hist.end(), 0 );
    uchar* p = hist.data();
    if( max > 0 ) {
        for( int i = 0; i < size; i++ ) {
            int R = 127 - ( 127 * r[i] ) / max;
            int G = 127 - ( 127 * g[i] / 2 ) / max;
            int B = 127 - ( 127 * b[i] ) / max;
            for( int j = 0; j < h; j++ ) {
                uchar* p0 = p + j * biteWidth + 3 * i;
                if( R < 127 ) {
                    if( j > R ) {
                        p0[0] = 0x50;
                    } else if( j == R ){
                        p0[0] = 0xFF;
                    }
                }
                if( G < 127 ) {
                    if( j > G ) {
                        p0[1] = 0x50;
                    } else if( j == G ){
                        p0[1] = 0xFF;
                    }
                }
                if( B < 127 ) {
                    if( j > B ) {
                        // Make blue on black brighter (otherwise it is too dim)
                        if( p0[0] == 0 && p0[1] == 0 ) {
                            p0[2] = 0x80;
                        } else {
                            p0[2] = 0x50;
                        }
                    } else if( j == B ){
                        // Make blue outline brighter
                        p0[0] = 0x80;
                        p0[1] = 0x80;
                        p0[2] = j < 127 ? 0xFF : 0;
                    }
                }
            }
        }
    }

    auto pixmap = Qt::CreatePixmap( p, size, h );

    // Add basic statistics (min/max etc)
    QPainter painter( &pixmap );
    QPen pen( QColor::fromRgb( 0x40, 0x80, 0x40 ) );
    painter.setPen( pen );
    QFont font( "Consolas" );
    font.setPointSizeF( 6.5 );
    painter.setFont( font );
    painter.drawText( size - 140, 10,
        QString( "max:%1 (%2%,%3)" ).arg(
            QString::number( maxValue ),
            QString::number( ( 100.0 * maxCount ) / ( width * height ), 'f', 1 ),
            collapseNumber( maxCount ) )
    );
    painter.drawText( size - 140, 20,
        QString( "min:%1 (%2%,%3)" ).arg(
            QString::number( minValue ),
            QString::number( ( 100.0 * minCount ) / ( width * height ), 'f', 1 ),
            collapseNumber( minCount ) )
    );
    painter.end();

    return pixmap;
}
