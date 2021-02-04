#include "imageview.h"

#include<QPainter>
#include "renderer.h"

#include <math.h>

#include <QDebug>

ImageView::ImageView(QWidget *parent) : QLabel(parent)
{

}

struct CChannelStat {
    int Median;
    int Sigma;
};

class CPixelStatistics {
public:
    CPixelStatistics( int numberOfChannels, int bitsPerChannel ) :
        channelSize( 2 << bitsPerChannel )
    {
        channels.resize( numberOfChannels );
        for( int i = 0; i < numberOfChannels; i++ ) {
            channels[i].resize( channelSize );
        }
    }

    std::vector<uint>& operator [] ( size_t index ) { return channels[index]; }
    const std::vector<uint>& operator [] ( size_t index ) const { return channels[index]; }

    CChannelStat stat( int channel, int count ) const
    {
        CChannelStat stat;
        stat.Median = p( channel, count / 2 );
        stat.Sigma = stat.Median - p( channel, count / 2 - count / 3 );
        int dR2 = p( channel, count / 2 + count / 3 ) - stat.Median;
        int mR = maxP( channel, stat.Median - stat.Sigma, stat.Median + dR2 );
        qDebug() << "Hist" << channel << stat.Median  << stat.Sigma << "/" << dR2 << mR;

        return stat;
    }

    size_t p( int channel, int target ) const;
    size_t maxP( int channel, int start, int end ) const;

private:
    const size_t channelSize;
    std::vector<std::vector<uint>> channels;
};

size_t CPixelStatistics::p( int channel, int target ) const
{
    auto& hist = channels[channel];
    uint sum = 0;
    for( int i = 0; i < channelSize; i++ ) {
        uint delta = target - sum;
        uint v = hist[i];
        if( delta <= v ) {
            if( i > 0 && delta <= v / 2 ) {
                return i - 1;
            } else {
                return i;
            }
        }
        sum += v;
    }
    return INT_MAX;
}

size_t CPixelStatistics::maxP( int channel, int start, int end ) const
{
    auto& hist = channels[channel];

    int max = INT_MIN;
    int maxPos = INT_MIN;
    for( int i = start; i <= end; i++ ) {
        int v = hist[i];
        if( v > max ) {
            max = v;
            maxPos = i;
        }
    }
    return maxPos;
}

void ImageView::render( int cx, int cy )
{
    Renderer renderer( image->RawPixels(), image->Width(), image->Height() );

    QPixmap pixmap = renderer.RenderRectHalfRes( cx, cy, 51, 51 );

    QPainter painter( &pixmap );
    QPen pen( QColor::fromRgb( 0xFF, 0, 0 ) );
    pen.setWidth( 1 );
    painter.setPen( pen );
    painter.drawEllipse( QPoint( cx, cy ), 5, 5 );

    setPixmap( pixmap );

    CPixelStatistics stats( 3, 12 );

    auto& histR = stats[0];
    auto& histG = stats[1];
    auto& histB = stats[2];

    const auto& size = histR.size();

    const ushort* raw = image->RawPixels();
    int width = image->Width();
    int height = image->Height();
    int W = 101;
    int H = 101;

    int count = 0;

    for( size_t y = cy - H / 2; y <= cy + H / 2; y++ ) {
        const ushort* srcLine = raw + 2 * width * y;
        for( size_t x = cx - W / 2; x <= cx + W / 2; x++ ) {
            const ushort* src = srcLine + 2 * x;
            ushort r = src[0];
            ushort g1 = src[1];
            ushort g2 = src[width];
            ushort b = src[width + 1];

            histR[r]++;
            histB[b]++;
            histG[g1]++;
            histG[g2]++;

            count++;
        }
    }

    const CChannelStat sR = stats.stat( 0, count );
    const CChannelStat sG = stats.stat( 1, count * 2);
    const CChannelStat sB = stats.stat( 2, count);

    size_t w = width / 2;
    size_t h = height / 2;
    size_t byteWidth = 3 * w;
    std::vector<uchar> pixels( byteWidth * h );
    uchar* rgb = pixels.data();
    for( size_t y = cy - H / 2; y <= cy + H / 2; y++ ) {
        const ushort* srcLine = raw + 2 * width * y;
        uchar* dstLine = rgb + byteWidth * y;
        for( size_t x = cx - W / 2; x <= cx + W / 2; x++ ) {
            const ushort* src = srcLine + 2 * x;
            uchar* dst = dstLine + 3 * x;

            uint r = src[0];
            uint g1 = src[1];
            uint g2 = src[width];
            uint b = src[width + 1];
            uint g = ( g1 + g2 ) / 2;

            const int k = 16;
            dst[0] = r <= ( sR.Median + k * sR.Sigma ) ? ( r < sR.Median ? 0 : ( 255 * ( r - sR.Median ) / k / sR.Sigma ) ) : 255;
            dst[1] = g <= ( sG.Median + k * sG.Sigma ) ? ( g < sG.Median ? 0 : ( 255 * ( g - sG.Median ) / k / sG.Sigma ) ) : 255;
            dst[2] = b <= ( sB.Median + k * sB.Sigma ) ? ( b < sB.Median ? 0 : ( 255 * ( b - sB.Median ) / k / sB.Sigma ) ) : 255;

            /*uint r = 0;
            uint g1 = 0;
            uint g2 = 0;
            uint b = 0;
            int count = 0;
            const int n = 5;
            for( int i = -n; i <= n; i++ )
            for( int j = -n; j <= n; j++ ) {
                int i0 = 2 * ( i * width + j );
                r += src[i0];
                g1 += src[i0 + 1];
                g2 += src[i0 + width];
                b += src[i0 + width + 1];
                count++;
            }
            r /= count;
            g1 /= count;
            g2 /= count;
            b /= count;
            uint g = ( g1 + g2 ) / 2;

            const int k = 4;
            const int m = 2 * n + 1;
            dst[0] = r < ( sR.Median + ( k * sR.Sigma ) / m ) ? 0 : 255;
            dst[1] = g < ( sG.Median + ( k * sG.Sigma ) / m ) ? 0 : 255;
            dst[2] = b < ( sB.Median + ( k * sB.Sigma ) / m ) ? 0 : 255;*/
        }
    }

    QImage image( rgb, w, h, QImage::Format_RGB888 );
    setPixmap( QPixmap::fromImage( image ) );

    size_t R = 0;
    for( R = 1; R < H / 2; R++ ) {
        size_t sqrR = R * R;
        size_t prevSqrR = ( R - 1 ) * ( R - 1);
        int sumR = 0;
        int sumSqrdR = 0;
        int count = 0;
        for( size_t y = cy - H / 2; y <= cy + H / 2; y++ ) {
            const ushort* srcLine = raw + 2 * width * y;
            for( size_t x = cx - W / 2; x <= cx + W / 2; x++ ) {
                size_t dx = x - cx;
                size_t dy = y - cy;
                size_t sqr = dx * dx + dy * dy;
                if( sqr <= sqrR && sqrR > prevSqrR ) {
                    const ushort* src = srcLine + 2 * x;
                    uint r = src[0];
                    uint g1 = src[1];
                    uint g2 = src[width];
                    uint b = src[width + 1];
                    uint g = ( g1 + g2 ) / 2;

                    sumR += r;
                    sumR -= sR.Median;
                    sumSqrdR += sR.Sigma * sR.Sigma;
                    count++;
                }
            }
        }
        if( sumR / count < 5 * sR.Sigma ) {
            R = R - 1;
            break;
        }
    }
    qDebug() << "R =" << R;

    double sumGdX = 0;
    double sumGdY = 0;
    double sumG = 0;
    for( size_t y = cy - R; y <= cy + R; y++ ) {
        int dy = y - cy;
        const ushort* srcLine = raw + 2 * width * y;
        uchar* dstLine = rgb + byteWidth * y;
        for( size_t x = cx - R; x <= cx + R; x++ ) {
            int dx = x - cx;
            if( dx * dx + dy * dy <= R * R ) {
                const ushort* src = srcLine + 2 * x;
                //uint r = src[0];
                uint g1 = src[1];
                uint g2 = src[width];
                //uint b = src[width + 1];
                //uint g = ( g1 + g2 ) / 2;

                double dg1 = g1 - sG.Median;
                if( dg1 > 0 ) {
                    sumGdX += dg1 * ( 2 * dx + 1 );
                    sumGdY += dg1 * 2 * dy;
                    sumG += dg1;
                }
                double dG2 = g2 - sG.Median;
                if( dG2 > 0 ) {
                    sumGdX += dG2 * 2 * dx;
                    sumGdY += dG2 * ( 2 * dy + 1 );
                    sumG += dG2;
                }
            }
        }
    }
    if( sumG > 0 ) {
        qDebug() << "dXg =" << - sumGdX / sumG;
        qDebug() << "dYg =" << - sumGdY / sumG;
    }

    double _cx = cx - sumGdX / sumG / 2;
    double _cy = cy - sumGdY / sumG / 2;
    
    cx = (int)round( _cx );
    cy = (int)round( _cy );

    double avgG = 0;
    uint countG = 0;
    for( size_t y = cy - 2 * R; y <= cy + 2 * R; y++ ) {
        int dy = y - cy;
        const ushort* srcLine = raw + 2 * width * y;
        uchar* dstLine = rgb + byteWidth * y;
        for( size_t x = cx - 2 * R; x <= cx + 2 * R; x++ ) {
            int dx = x - cx;
            int rr = dx * dx + dy * dy;
            if( rr > R * R && rr <= 4 * R * R ) {
                const ushort* src = srcLine + 2 * x;
                uint r = src[0];
                uint g1 = src[1];
                uint g2 = src[width];
                uint b = src[width + 1];
                uint g = ( g1 + g2 ) / 2;
                //if( g < sG.Median + 5 * sG.Sigma ) {
                    avgG += g;
                    countG++;
                //}
            }
        }
    }
    avgG = avgG / countG;

    double sumGD = 0;
    double _sumG = 0;
    int maxG = 0;
    int maxCount = 0;
    int minG = INT_MAX;
    for( int y = cy - R; y <= cy + R; y++ ) {
        double dy = y - _cy;
        const ushort* srcLine = raw + 2 * width * y;
        uchar* dstLine = rgb + byteWidth * y;
        for( int x = cx - R; x <= cx + R; x++ ) {
            double dx = x - _cx;
            if( dx * dx + dy * dy <= R * R ) {
                const ushort* src = srcLine + 2 * x;
                uint r = src[0];
                uint g1 = src[1];
                uint g2 = src[width];
                uint b = src[width + 1];
                uint g = ( g1 + g2 ) / 2;

                uchar* dst = dstLine + 3 * x;

                double dG = g1 - avgG;
                if( dG > 3 * sG.Sigma ) {
                    double d = sqrt( ( 2 * dx + 1 ) * ( 2 * dx + 1 ) + 4 * dy * dy );
                    sumGD += dG * d;
                    _sumG += dG;

                    dst[0] = 0;
                    dst[1] = dG > 0 ? 255 : 0;
                    dst[2] = 0;
                }
                dG = g2 - avgG;
                if( dG > 3 * sG.Sigma ) {
                    double d = sqrt( 4 * dx * dx + ( 2 * dy + 1 ) * ( 2 * dy + 1 ) );
                    sumGD += dG * d;
                    _sumG += dG;

                    dst[0] = 0;
                    dst[1] = dG > 0 ? 255 : 0;
                    dst[2] = 0;
                }
                if( g >= maxG ) {
                    if( g > maxG ) {
                        maxG = g;
                        maxCount = 0;
                    }
                    maxCount++;
                }
                if( g < minG ) {
                    minG = g;
                }
            }
        }
    }

    if( sumG > 0 ) {
        qDebug() << "Gmin/max =" << minG << maxG;
        qDebug() << "maxCount =" << maxCount;
        qDebug() << "HFD =" << ( 1.0 * sumGD ) / _sumG;
        qDebug() << "Flux =" << _sumG;
    }

    QImage image2( rgb, w, h, QImage::Format_RGB888 );
    setPixmap( QPixmap::fromImage( image2 ) );
    
    R *= 2;
    std::vector<double> pD( R );
    std::vector<uint> counts( R );
    for( size_t y = cy - R; y <= cy + R; y++ ) {
        int dy = y - cy;
        const ushort* srcLine = raw + 2 * width * y;
        for( size_t x = cx - R; x <= cx + R; x++ ) {
            int dx = x - cx;
            int rr = (int)std::round( std::sqrt( dx * dx + dy * dy ) );
            if( rr < R ) {
                const ushort* src = srcLine + 2 * x;
                uint r = src[0];
                uint g1 = src[1];
                uint g2 = src[width];
                uint b = src[width + 1];
                uint g = ( g1 + g2 ) / 2;
                
                pD[rr] += g;
                counts[rr]++;
            }
        }
    }

    double maxV = 0;
    for( int i = 0; i < R; i++ ) {
        double v = pD[i] / counts[i] - avgG;
        if( v > maxV ) {
            maxV = v;
        }
        pD[i] = v;
    }

    w = R;
    h = 128;
    int fh = h + 5;
    byteWidth = R * 3;

    std::vector<uchar> hist( fh * byteWidth );
    std::fill( hist.begin(), hist.end(), 0 );
    uchar* p = hist.data();
    for( int i = 0; i < w; i++ ) {
        double _v = h - ( h * pD[i] ) / maxV;
        int v = (int)( round( _v ) );
        for( int j = 0; j < fh; j++ ) {
            if( j >= v ) {
                 uchar* p0 = p + j * byteWidth + 3 * i;
                 if( v < h && j < h ) {
                     p0[0] = 0xFF;
                     p0[1] = 0xFF;
                     p0[2] = 0xFF;
                 } else if ( v > h && j > h ) {
                     p0[0] = 0xFF;
                     p0[1] = 0;
                     p0[2] = 0;
                 } else if( v == h && j == h ) {
                     if( _v <= v ) {
                         char c = (char)( 256 * (_v - v ) );
                         p0[0] = c;
                         p0[1] = c;
                         p0[2] = c;
                     } else {
                         char c = (char)( 256 * (v - _v ) );
                         p0[0] = c;
                         p0[1] = 0;
                         p0[2] = 0;
                     }
                 }
            }
        }
    }

    QImage image3( p, w, fh, byteWidth, QImage::Format_RGB888 );
    setPixmap( QPixmap::fromImage( image3 ) );
}
