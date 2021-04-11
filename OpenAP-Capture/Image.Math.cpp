// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#include "Image.Math.h"

#include <Image.Debayer.HQLinear.h>

#include <Math.Geometry.h>
#include <Math.LinearAlgebra.h>

#include <QDebug>

#include <math.h>
#include <stack>
#include <algorithm>

CPixelStatistics::CPixelStatistics( int numberOfChannels, int bitsPerChannel ) :
    channelSize( 2 << bitsPerChannel )
{
    channels.resize( numberOfChannels );
    for( int i = 0; i < numberOfChannels; i++ ) {
        channels[i].resize( channelSize );
    }
}

CChannelStat CPixelStatistics::stat( int channel, int n ) const
{
    CChannelStat stat;
    int count = n * this->count;
    stat.Median = p( channel, count / 2 );
    stat.Sigma = stat.Median - p( channel, count / 2 - count / 3 );
    int dR2 = p( channel, count / 2 + count / 3 ) - stat.Median;
    int mR = maxP( channel, stat.Median - stat.Sigma, stat.Median + dR2 );
    return stat;
}

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

CRawU16::CRawU16( const CRawU16Image* image ) :
    CRawU16( image->RawPixels(), image->Width(), image->Height(), image->BitDepth() )
{
}

CRawU16::CRawU16( const ushort* _raw, int _width, int _height, int _bitDepth ) :
    raw( _raw ), width( _width ), height( _height ), bitDepth( _bitDepth )
{

}

std::shared_ptr<CRgbU16Image> CRawU16::DebayerRect( int x, int y, int w, int h ) const
{
    auto result = std::make_shared<CRgbU16Image>( w, h );
    CDebayer_RawU16_HQLinear debayer( raw, width, height, bitDepth );
    debayer.ToRgbU16( result->RgbPixels(), result->Stride(), x, y, w, h );
    return result;
}

std::shared_ptr<CGrayU16Image> CRawU16::GrayU16( int x, int y, int width, int height ) const
{
    return ToGrayU16( DebayerRect( x, y, width, height ).get() );
}

CPixelStatistics CRawU16::CalculateStatistics( int x0, int y0, int W, int H ) const
{
    CPixelStatistics stats( 3, bitDepth );

    auto& histR = stats[0];
    auto& histG = stats[1];
    auto& histB = stats[2];

    int count = 0;

    int cx = x0 + W / 2;
    int cy = y0 + H / 2;
    x0 = cx - W;
    y0 = cy - H;
    x0 += x0 % 2;
    y0 += y0 % 2;

    for( size_t y = 0; y < H; y++ ) {
        int Y = 2 * y + y0;
        if( Y < 0 || Y >= height ) {
            continue;
        }
        const ushort* srcLine = raw + width * Y;
        for( size_t x = 0; x < W; x++ ) {
            int X = 2 * x + x0;
            if( X < 0 || X >= width ) {
                continue;
            }
            const ushort* src = srcLine + X;
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

    stats.setCount( count );

    return stats;
}

std::shared_ptr<CRgbImage> CRawU16::Stretch( int x0, int y0, int W, int H ) const
{
    CPixelStatistics stats = CalculateStatistics( x0, y0, W, H );

    const int maxValue = ~(~0u << bitDepth) - 1;

    const CChannelStat sR = stats.stat( 0 );
    const CChannelStat sG = stats.stat( 1, 2);
    const CChannelStat sB = stats.stat( 2 );

    auto rgb16 = DebayerRect( x0, y0, W, H );

    auto result = std::make_shared<CRgbImage>( W, H );
    for( size_t y = 0; y < H; y++ ) {
        const ushort* srcLine = rgb16->ScanLine( y );
        uchar* dstLine = result->ScanLine( y );
        for( size_t x = 0; x < W; x++ ) {
            const ushort* src = srcLine + 3 * x;
            uchar* dst = dstLine + 3 * x;

            uint r = src[0];
            uint g = src[1];
            uint b = src[2];

            if( r >= maxValue || g >= maxValue || b >= maxValue ) {
                dst[0] = 0xFF;
                dst[1] = 0x00;
                dst[2] = 0x80;
            } else {
                const int k = 12;
                dst[0] = r <= ( sR.Median + k * sR.Sigma ) ? ( r < sR.Median ? 0 : ( 255 * ( r - sR.Median ) / k / sR.Sigma ) ) : 255;
                dst[1] = g <= ( sG.Median + k * sG.Sigma ) ? ( g < sG.Median ? 0 : ( 255 * ( g - sG.Median ) / k / sG.Sigma ) ) : 255;
                dst[2] = b <= ( sB.Median + k * sB.Sigma ) ? ( b < sB.Median ? 0 : ( 255 * ( b - sB.Median ) / k / sB.Sigma ) ) : 255;
            }
        }
    }

    return result;
}

std::shared_ptr<CRgbImage> CRawU16::StretchHalfRes( int x0, int y0, int W, int H ) const
{
    CPixelStatistics stats = CalculateStatistics( x0, y0, W, H );

    const int maxValue = ~(~0u << bitDepth) - 1;

    const CChannelStat sR = stats.stat( 0 );
    const CChannelStat sG = stats.stat( 1, 2);
    const CChannelStat sB = stats.stat( 2 );

    W /= 2;
    H /= 2;

    auto result = std::make_shared<CRgbImage>( W, H );
    for( size_t y = 0; y < H; y++ ) {
            int Y = 2 * y + y0;
            if( Y < 0 || Y >= height ) {
                continue;
            }
            const ushort* srcLine = raw + width * Y;
            uchar* dstLine = result->RgbPixels() + result->ByteWidth() * y;
            for( size_t x = 0; x < W; x++ ) {
                int X = 2 * x + x0;
                if( X < 0 || X >= width ) {
                    continue;
                }
                const ushort* src = srcLine + X;
                uchar* dst = dstLine + 3 * x;

                uint r = src[0];
                uint g1 = src[1];
                uint g2 = src[width];
                uint b = src[width + 1];
                uint g = ( g1 + g2 ) / 2;

                if( r >= maxValue || g1 >= maxValue || g2 >= maxValue || b >= maxValue ) {
                    dst[0] = 0xFF;
                    dst[1] = 0x00;
                    dst[2] = 0x80;
                } else {
                    const int k = 12;
                    dst[0] = r <= ( sR.Median + k * sR.Sigma ) ? ( r < sR.Median ? 0 : ( 255 * ( r - sR.Median ) / k / sR.Sigma ) ) : 255;
                    dst[1] = g <= ( sG.Median + k * sG.Sigma ) ? ( g < sG.Median ? 0 : ( 255 * ( g - sG.Median ) / k / sG.Sigma ) ) : 255;
                    dst[2] = b <= ( sB.Median + k * sB.Sigma ) ? ( b < sB.Median ? 0 : ( 255 * ( b - sB.Median ) / k / sB.Sigma ) ) : 255;
                }
            }
        }

    return result;
}

std::shared_ptr<CRgbImage> CRawU16::StretchQ( int x, int y, int w, int h ) const
{
    auto rgbImage = StretchHalfRes( x, y, w, h );
    int byteWidth = rgbImage->ByteWidth();
    w /= 2;
    h /= 2;

    auto result = std::make_shared<CRgbImage>( w / 2, h / 2 );
    for( int i = 0; i < h / 2; i++ ) {
        const uchar* ptr = rgbImage->ScanLine( 2 * i );
        uchar* ptr2 = result->ScanLine( i );
        for( int j = 0; j < w / 2; j++ ) {
            const uchar* src0 = ptr + 2 * 3 * j;
            const uchar* src1 = src0 + 3;
            const uchar* src2 = src0 + byteWidth;
            const uchar* src3 = src0 + byteWidth + 3;

            uchar* dst = ptr2 + 3 * j;
            dst[0] = ( src0[0] + src1[0] + src2[0] + src3[0] ) / 4;
            dst[1] = ( src0[1] + src1[1] + src2[1] + src3[1] ) / 4;
            dst[2] = ( src0[2] + src1[2] + src2[2] + src3[2] ) / 4;
        }
    }
    return result;
}

std::shared_ptr<CGrayU16Image> CRawU16::ToGrayU16( const CRgbU16Image* rgb16 )
{
    int width = rgb16->Width();
    int height = rgb16->Height();
    auto result = std::make_shared<CGrayU16Image>( width, height );
    for( int y = 0; y < height; y++ ) {
        const ushort* srcLine = rgb16->ScanLine( y );
        ushort* dstLine = result->ScanLine( y );
        for( int x = 0; x < width; x++ ) {
            const ushort* src = srcLine + 3 * x;
            uint r = src[0];
            uint g = src[1];
            uint b = src[2];

            dstLine[x] = r + g + b;
        }
    }

    return result;
}

std::shared_ptr<CGrayImage> CRawU16::ToGray( const CGrayU16Image* rgb16 )
{
    int width = rgb16->Width();
    int height = rgb16->Height();
    auto result = std::make_shared<CGrayImage>( width, height );
    for( int y = 0; y < height; y++ ) {
        const ushort* srcLine = rgb16->ScanLine( y );
        uchar* dstLine = result->ScanLine( y );
        for( int x = 0; x < width; x++ ) {
            const ushort v = srcLine[x] >> 3;
            dstLine[x] = v < 255 ? v : 255;
        }
    }

    return result;
}

void CRawU16::GradientAscentToLocalMaximum( int& x, int& y, int size )
{
    auto gray = GrayU16( x - size / 2, y - size / 2, size, size );
    int _x = size / 2;
    int _y = size / 2;
    CRawU16::GradientAscentToLocalMaximum( gray.get(), _x, _y, 25 );
    CRawU16::GradientAscentToLocalMaximum( gray.get(), _x, _y, 3 );
    x += _x - size / 2;
    y += _y - size / 2;
}

void CRawU16::GradientAscentToLocalMaximum( const CGrayU16Image* image, int& x, int& y, int window )
{
    int halfW = window / 2;
    int stride = image->Stride();
    while( true ) {
        const ushort* p = image->ScanLine( y ) + x;
        int maxv = 0;
        int maxPos = 0;
        int pos = 0;
        // For each pixel around the current pixel
        for( int i = -halfW; i <= halfW; i+= halfW ) {
            for( int j = -halfW; j <= halfW; j+= halfW ) {
                const ushort* s = p + i * stride + j;
                int v = 0;
                // Sum values of pixels in 3x3 area
                for( int n = -halfW; n <= halfW; n++ ) {
                    for( int m = -halfW; m <= halfW; m++ ) {
                        v += s[n * stride + m];
                    }
                }
                if( v > maxv ) {
                    maxv = v;
                    maxPos = pos;
                }
                pos++;
            }
        }
        // Move to the found maximum
        switch( maxPos ) {
            case 0: x--; y--; break;
            case 1: y--; break;
            case 2: x++; y--; break;
            case 3: x--; break;
            case 4: return;
            case 5: x++; break;
            case 6: x--; y++; break;
            case 7: y++; break;
            case 8: x++; y++; break;
        }
    }
}

CPixelStatistics CRawU16::CalculateStatistics( const CGrayU16Image* image )
{
    CPixelStatistics stats( 1, 16 );
    auto& hist = stats[0];

    int width = image->Width();
    int height = image->Height();

    int count = 0;
    for( int y = 0; y < height; y++ ) {
        const ushort* srcLine = image->ScanLine( y );
        for( int x = 0; x < width; x++ ) {
            const ushort* src = srcLine + x;
            ushort v = src[0];
            hist[v]++;
            count++;
        }
    }
    stats.setCount( count );

    return stats;
}

static std::tuple<int, int, int, int, int, unsigned> starMaskWithMax( const CGrayU16Image* image, int x, int y, int t, int width, int height, std::shared_ptr<CGrayImage> result, int value )
{
    int vmax = 0;
    int xmax = 0;
    int ymax = 0;
    int cmax = 0;
    int count = 0;
    unsigned sum = 0;
    std::stack<std::tuple<int, int>> s;
    s.emplace( x, y );
    while( not s.empty() ) {
        int x, y;
        std::tie( x, y ) = s.top();
        s.pop();
        if( x >= 0 && x < width && y > 0 && y < height ) {
            auto& dst = result->At( x, y );
            int v = image->At( x, y );
            if( dst < value && v >= t ) {
                dst = value;
                count++;
                sum += v;
                if( v >= vmax ) {
                    if( v > vmax ) {
                        vmax = v;
                        xmax = x;
                        ymax = y;
                        cmax = 1;
                    } else {
                        cmax++;
                    }
                }
                s.emplace( x - 1, y );
                s.emplace( x + 1, y );
                s.emplace( x - 1, y - 1 );
                s.emplace( x, y - 1);
                s.emplace( x + 1, y - 1);
                s.emplace( x - 1, y + 1 );
                s.emplace( x, y + 1);
                s.emplace( x + 1, y + 1);
            }
        }
    }
    return std::make_tuple( vmax, xmax, ymax, cmax, count, sum );
}

static void erase( int x, int y, int width, int height, std::shared_ptr<CGrayImage> result )
{
    std::stack<std::tuple<int, int>> s;
    s.emplace( x, y );
    while( not s.empty() ) {
        int x, y;
        std::tie( x, y ) = s.top();
        s.pop();
        if( x >= 0 && x < width && y > 0 && y < height ) {
            auto& dst = result->At( x, y );
            if( dst != 0  ) {
                dst = 0;
                s.emplace( x - 1, y );
                s.emplace( x + 1, y );
                s.emplace( x - 1, y - 1 );
                s.emplace( x, y - 1);
                s.emplace( x + 1, y - 1);
                s.emplace( x - 1, y + 1 );
                s.emplace( x, y + 1);
                s.emplace( x + 1, y + 1);
            }
        }
    }
}

static void enumerate(const CGrayU16Image* image, int x, int y, int width, int height, int th, std::shared_ptr<CGrayImage> mask, std::function<void(int, int, ushort)> f )
{
    std::stack<std::tuple<int, int, int>> undo;
    std::stack<std::tuple<int, int>> s;
    s.emplace( x, y );
    while( not s.empty() ) {
        int x, y;
        std::tie( x, y ) = s.top();
        s.pop();
        if( x >= 0 && x < width && y > 0 && y < height ) {
            auto& dst = mask->At( x, y );
            if( dst >= th  ) {
                f( x, y, image->At( x, y ) );
                undo.emplace( x, y, dst );
                dst = 0;
                s.emplace( x - 1, y );
                s.emplace( x + 1, y );
                s.emplace( x - 1, y - 1 );
                s.emplace( x, y - 1);
                s.emplace( x + 1, y - 1);
                s.emplace( x - 1, y + 1 );
                s.emplace( x, y + 1);
                s.emplace( x + 1, y + 1);
            }
        }
    }
    while( not undo.empty() ) {
        int x, y, v;
        std::tie( x, y, v ) = undo.top();
        undo.pop();
        mask->At( x, y ) = v;
    }
}

static std::shared_ptr<CGrayImage> starMask( const CGrayU16Image* image, int x, int y, int t, int width, int height, std::shared_ptr<CGrayImage> result, int value )
{
    std::stack<std::tuple<int, int>> s;
    s.emplace( x, y );
    while( not s.empty() ) {
        int x, y;
        std::tie( x, y ) = s.top();
        s.pop();
        if( x >= 0 && x < width && y > 0 && y < height ) {
            auto& dst = result->At( x, y );
            if( dst < value && image->At( x, y ) >= t ) {
                dst = value;
                s.emplace( x - 1, y );
                s.emplace( x + 1, y );
                s.emplace( x - 1, y - 1 );
                s.emplace( x, y - 1);
                s.emplace( x + 1, y - 1);
                s.emplace( x - 1, y + 1 );
                s.emplace( x, y + 1);
                s.emplace( x + 1, y + 1);
            }
        }
    }
    return result;
}

static std::shared_ptr<CGrayImage> starMask( const CGrayU16Image* image, int x, int y, int t, int value = 128 )
{
    int width = image->Width();
    int height = image->Height();
    std::shared_ptr<CGrayImage> result = std::make_shared<CGrayImage>( width, height );
    return starMask( image, x, y, t, width, height, result, value );
}

static std::shared_ptr<CGrayImage> starMask( const CGrayU16Image* image, int x, int y, int t, int value, std::shared_ptr<CGrayImage>& result )
{
    int width = image->Width();
    int height = image->Height();
    return starMask( image, x, y, t, width, height, result, value );
}

class CalculateXY {
public:
    CalculateXY( int _x0, int _y0, int _halfMax ) : x0( _x0 ), y0( _y0 ), halfMax( _halfMax ){}
    void operator()( int x, int y, ushort v )
    {
        double value = v - halfMax;
        if( value > 0 ) {
            sumVX += value * ( x - x0 );
            sumVY += value * ( y - y0 );
            sumV += value;
        }
    }

    double dx() { return sumVX / sumV; }
    double dy() { return sumVY / sumV; }

private:
    int x0;
    int y0;
    int halfMax;
    double sumVX = 0;
    double sumVY = 0;
    double sumV = 0;
};

struct StarXY {
    double NormalizedBrightness;
    double X;
    double Y;
    std::vector<int> Neighbours;
    std::vector<double> Distances;
    std::vector<double> Angles;
    int MatchedIndex = -1;
};

static std::vector<StarXY> normalizedBrightStarsXY( const DetectionResults& r0, size_t num )
{
    auto d0 = r0.DetectionRegions;
    auto sortByBrightness = [&]( std::shared_ptr<DetectionRegion> r0, std::shared_ptr<DetectionRegion> r1 ) { return r0->FluxAtHalfDetectionThreshold > r1->FluxAtHalfDetectionThreshold; };
    std::sort( d0.begin(), d0.end(), sortByBrightness );
    auto size = std::min( d0.size(), num );
    auto b0 =  d0.front()->FluxAtHalfDetectionThreshold;

    std::vector<StarXY> result;
    result.resize( size );
    for( size_t i = 0; i < size; i++ ) {
        auto dr0 = *d0[i];
        result[i].NormalizedBrightness = 1.0 * dr0.FluxAtHalfDetectionThreshold / b0;
        CalculateXY c0( dr0.Xmax, dr0.Ymax, dr0.Background + ( dr0.Vmax - dr0.Background ) / 2 );
        enumerate( r0.Image.get(), dr0.Xmax, dr0.Ymax, r0.Image->Width(), r0.Image->Height(), 128, r0.Mask, std::ref( c0 ) );
        result[i].X = dr0.Xmax + c0.dx();
        result[i].Y = dr0.Ymax + c0.dy();
    }
    return result;
}

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    const std::vector<T>& vec,
    const Compare& compare)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
    return p;
}

template <typename T>
std::vector<T> apply_permutation(
    const std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](std::size_t i){ return vec[i]; });
    return sorted_vec;
}

static void calculateDistances( std::vector<StarXY>& stars )
{
    for( size_t i = 0; i < stars.size(); i++ ) {
        auto& d = stars[i].Distances;
        std::vector<double> DX;
        std::vector<double> DY;
        auto& neighbours = stars[i].Neighbours;
        for( size_t j = 0; j < i; j++ ) {
            double dX = stars[j].X - stars[i].X;
            double dY = stars[j].Y - stars[i].Y;
            double L = sqrt( dX * dX + dY * dY );
            d.push_back( L );
            DX.push_back( dX );
            DY.push_back( dY );
            neighbours.push_back( j );
        }
        if( d.size() < 7 ) {
            for( size_t j = i + 1; j < stars.size() && d.size() < 7; j++ ) {
                double dX = stars[j].X - stars[i].X;
                double dY = stars[j].Y - stars[i].Y;
                double L = sqrt( dX * dX + dY * dY );
                d.push_back( L );
                DX.push_back( dX );
                DY.push_back( dY );
                neighbours.push_back( j );
            }
        }
        if( d.size() > 7 ) {
            auto p = sort_permutation( d, std::less<double>() );
            d = apply_permutation( d, p );
            DX = apply_permutation( DX, p );
            DY = apply_permutation( DY, p );
            neighbours = apply_permutation( neighbours, p );
            d.resize( 7 );
            DX.resize( 7 );
            DY.resize( 7 );
            neighbours.resize( 7 );

            std::vector<double> brightness;
            for( size_t j = 0; j < neighbours.size(); j++ ) {
                brightness.push_back( stars[neighbours[j]].NormalizedBrightness );
            }
            auto p1 = sort_permutation( brightness, std::greater<double>() );
            d = apply_permutation( d, p1 );
            DX = apply_permutation( DX, p1 );
            DY = apply_permutation( DY, p1 );
            neighbours = apply_permutation( neighbours, p1 );
        }
        auto& a = stars[i].Angles;
        for( size_t j = 0; j < d.size(); j++ ) {
            a.push_back( atan2( DX[j], DY[j] ) );
        }
    }
}

static bool match( const StarXY& star0, const StarXY& star1 )
{
    const auto& n0 = star0.Neighbours;
    const auto& n1 = star1.Neighbours;
    const auto& d0 = star0.Distances;
    const auto& d1 = star1.Distances;
    const auto& a0 = star0.Angles;
    const auto& a1 = star1.Angles;
    int count = 0;
    int fCount = 0;
    for( size_t n = 0; n < std::min( n0.size(), n1.size() ); n++ ) {
        double L0 = d0[n];
        double A0 = a0[n] - a0[0];
        double L1 = d1[n];
        double A1 = a1[n] - a1[0];
        if( fabs( L0 - L1 ) < 2 && fabs( A0 - A1 ) < 2 / L0 ) {
            count++;
        } else {
            fCount++;
        }
    }
    if( fCount == 0 || count > 2 ) {
        return true;
    }
    count = 0;
    for( size_t n = 0; n < std::min( n0.size(), n1.size() ); n++ ) {
        double L0 = d0[n];
        double A0 = a0[n] - a0[1];
        double L1 = d1[n];
        double A1 = a1[n] - a1[1];
        if( fabs( L0 - L1 ) < 2 && fabs( A0 - A1 ) < 2 / L0 ) {
            count++;
        }
    }
    return count > 2;
}

std::vector<std::pair<double, double>> align( std::vector<StarXY>& stars0, std::vector<StarXY>& stars1 )
{
    calculateDistances( stars0 );
    calculateDistances( stars1 );

    std::vector<std::pair<double, double>> result;
    for( size_t i = 0; i < stars0.size(); i++ ) {
        const auto& s0 = stars0[i];
        bool matched = false;
        for( size_t j = 0; j < stars1.size(); j++ ) {
            auto& s1 = stars1[j];
            if( match( s0, s1 ) ) {
                qDebug() << "Matched" << i << j << s0.X << s0.Y << s1.X << s1.Y;
                assert( !matched );
                matched = true;
                result.push_back( std::make_pair( s1.X, s1.Y ) );
                s1.MatchedIndex = i;
            }
        }
        if( !matched ) {
            result.push_back( std::make_pair( 0, 0 ) );
        }
    }
    return result;
}

void CFocusingHelper::AddFrame( const CRawU16Image* currentImage, int imageSize, int focuserPos )
{
    double prevCX = CX;
    double prevCY = CY;

    // Lock on the star (center on local maximum)
    CRawU16 rawU16( currentImage );
    rawU16.GradientAscentToLocalMaximum( cx, cy, imageSize );
    auto image = rawU16.GrayU16( cx - imageSize / 2, cy - imageSize / 2, imageSize, imageSize );

    int height = image->Height();
    int width = image->Width();

    const auto s = CRawU16::CalculateStatistics( image.get() ).stat( 0 );
    qDebug() << "Median:" << s.Median << "Sigma:" << s.Sigma;

    int cVal = image->At( imageSize / 2, imageSize / 2 );

    int maxVal = 0;
    int halfMax = s.Median + ( cVal - s.Median ) / 2;
    int count = 0;
    for( int y = imageSize / 2 - 10; y <= imageSize / 2 + 10; y++ ) {
        const ushort* srcLine = image->ScanLine( y );
        for( int x = imageSize / 2 - 10; x < imageSize / 2 + 10; x++ ) {
            const ushort* src = srcLine + x;
            ushort v = src[0];
            if( v >= halfMax ) {
                count++;
            }
            if( v > maxVal ) {
                maxVal = v;
            }
        }
    }

    qDebug() << "Value at center:" << cVal;
    qDebug() << "Max val:" << maxVal;
    qDebug() << "Count at half max:" << count;
    qDebug() << "FWHM:" << 2 * sqrt ( count / M_PI );


    auto mask = starMask( image.get(), imageSize / 2, imageSize / 2, s.Median + s.Sigma, 16 );
    starMask( image.get(), imageSize / 2, imageSize / 2, s.Median + ( cVal - s.Median ) / 10, 32, mask );
    starMask( image.get(), imageSize / 2, imageSize / 2, s.Median + ( cVal - s.Median ) / 2, 128, mask );
    starMask( image.get(), imageSize / 2, imageSize / 2, s.Median + ( 9 * ( cVal - s.Median ) ) / 10, 192, mask );
    starMask( image.get(), imageSize / 2, imageSize / 2, s.Median + ( 19 * ( cVal - s.Median ) ) / 20, 255, mask );

    double sumVX = 0;
    double sumVY = 0;
    double sumV = 0;
    for( int y = 0; y < height; y++ ) {
        const ushort* src = image->ScanLine( y );
        const uchar* msk = mask->ScanLine( y );
        for( int x = 0; x < width; x++ ) {
            if( msk[x] > 0 ) {
                double value = src[x] - halfMax;
                if( value > 0 ) {
                    sumVX += value * ( x - imageSize / 2 );
                    sumVY += value * ( y - imageSize / 2 );
                    sumV += value;
                }
            }
        }
    }
    double dX = sumVX / sumV;
    double dY = sumVY / sumV;
    qDebug() << "dX:" << QString::number( dX, 'f', 2 ) << "dY:" << QString::number( dY, 'f', 2 );

    double meanSky = INT_MAX;
    if( R == 0 ) {
        const int dR = 10;
        int bestR = 0;
        int samples = 0;
        while( true ) {
            double sumV = 0;
            int count = 0;
            for( int y = 0; y < height; y++ ) {
                const ushort* src = image->ScanLine( y );
                for( int x = 0; x < width; x++ ) {
                    int dx = x - imageSize / 2;
                    int dy = y - imageSize / 2;
                    int rr = dx * dx + dy * dy;
                    if( rr >= R * R && rr < ( R + dR ) * ( R + dR ) ) {
                        sumV += src[x];
                        count++;
                    }
                }
            }
            double current = sumV / count;
            //qDebug() << "R:" << R << "Mean Sky:" << current;
            bool foundBoundary = false;
            if( current < meanSky || R == 0 ) {
                meanSky = current;
                bestR = R;
            } else {
                foundBoundary = true;
                //qDebug() << "----------";
            }
            if( foundBoundary || samples > 0 ) {
                if( ++samples == 15 ) {
                    break;
                }
            }
            R += 3;
        }
        R = bestR;
        R_out = R + dR;
    } else {
        double sumV = 0;
        int count = 0;
        for( int y = 0; y < height; y++ ) {
            const ushort* src = image->ScanLine( y );
            for( int x = 0; x < width; x++ ) {
                int dx = x - imageSize / 2;
                int dy = y - imageSize / 2;
                int rr = dx * dx + dy * dy;
                if( rr >= R * R && rr < ( R + 10 ) * ( R + 10 ) ) {
                    sumV += src[x];
                    count++;
                }
            }
        }
        meanSky = sumV / count;
    }
    qDebug() << "R:" << R << "Mean Sky:" << meanSky;

    sumV = 0;
    double sumVR = 0;
    for( int y = 0; y < height; y++ ) {
        const ushort* src = image->ScanLine( y );
        const uchar* msk = mask->ScanLine( y );
        for( int x = 0; x < width; x++ ) {
            if( msk[x] >= 32 ) {
                double value = src[x] - meanSky; // (?) ( s.Median + ( cVal - s.Median ) / 10 );
                for( int n = 0; n < 3; n++ ) {
                    for( int m = 0; m < 3; m++ ) {
                        double dx = x - ( imageSize / 2 + dX ) + n / 3.0;
                        double dy = y - ( imageSize / 2 + dY ) + m / 3.0;
                        double r = sqrt( dx * dx + dy * dy );
                        if( r < R ) {
                            sumVR += value * r;
                            sumV += value;
                        }
                    }
                }
            }
        }
    }
    HFD = 2 * sumVR / sumV;
    qDebug() << "HFD:" << QString::number( HFD, 'f', 2 );

    double halfMax2 = meanSky + ( maxVal - meanSky ) / 2;
    count = 0;
    for( int y = imageSize / 2 - 10; y <= imageSize / 2 + 10; y++ ) {
        const ushort* src = image->ScanLine( y );
        for( int x = imageSize / 2 - 10; x < imageSize / 2 + 10; x++ ) {
            if( src[x] - 0.5 >= halfMax2 ) {
                count++;
            }
        }
    }

    qDebug() << "(2) Count at half max:" << count;
    qDebug() << "(2) FWHM:" << 2 * sqrt ( count / M_PI );
    qDebug() << "(2) Max:" << maxVal;

    Mask = mask;

    CX = cx + dX;
    CY = cy + dY;

    auto it = focuserStats.find( focuserPos );
    if( it == focuserStats.end() ) {
        currentSeries = std::make_shared<Data>();
        focuserStats[focuserPos] = currentSeries;
        focuserPositions.push_back( focuserPos );
        std::sort( focuserPositions.begin(), focuserPositions.end() );
    } else {
        currentSeries = it->second;
    }

    currentSeries->HFD.push_back( HFD );
    currentSeries->FWHM.push_back( 2 * sqrt ( count / M_PI ) );
    currentSeries->Max.push_back( maxVal );
    currentSeries->CX.push_back( CX );
    currentSeries->CY.push_back( CY );

    if( isGlobalPolarAlign ) {
        currentSeries->theDetectionResults.push_back( rawU16.DetectStars( 0, 0, currentImage->Width(), currentImage->Height() ) );
    } else {
        currentSeries->theDetectionResults.push_back( DetectionResults() );
    }

    double sumdCX = 0;
    double sumdCY = 0;
    double sumdCXdCX = 0;
    double sumdCYdCY = 0;
    dCX = CX - prevCX;
    dCY = CY - prevCY;

    int minSize = currentSeries->CX.size();
    for( auto helper : extra ) {
        if( helper->currentSeries ) {
            int size = helper->currentSeries->CX.size();
            if( size < minSize ) {
                minSize = size;
            }
        }
    }

    if( extra.size() > 0 ) {
        sumdCX += dCX;
        sumdCXdCX += dCX * dCX;
        sumdCY += dCY;
        sumdCYdCY += dCY * dCY;
        int countC = 1;

        CFocusingHelper* prev = this;
        double d = 0; 
        for( auto helper : extra ) {
            double prevCX = helper->CX;
            double prevCY = helper->CY;
            helper->AddFrame( currentImage, imageSize, focuserPos );

            double dCX = helper->CX - prevCX;
            double dCY = helper->CY - prevCY;
            sumdCX += dCX;
            sumdCXdCX += dCX * dCX;
            sumdCY += dCY;
            sumdCYdCY += dCY * dCY;
            countC++;

            double dx = prev->CX - helper->CX;
            double dy = prev->CY - helper->CY;
            d += sqrt( dx * dx + dy * dy );
            prev = helper.get();
        }      
        sumD += d;
        sumDD += d * d;
        countL++;
        L = sumD / countL;
        sigmaL = sqrt( sumDD / countL - ( sumD * sumD ) / countL / countL );

        sigmadCX = sqrt( sumdCXdCX / countC - ( sumdCX * sumdCX ) / countC / countC );
        sigmadCY = sqrt( sumdCYdCY / countC - ( sumdCY * sumdCY ) / countC / countC );

        dCX = sumdCX / countC;
        dCY = sumdCY / countC;

        if( extra.size() > 1 ) {
            if( minSize >= 2 ) {
                std::vector<double> x1;
                std::vector<double> y1;
                x1.push_back( currentSeries->CX.cend()[-minSize] );
                y1.push_back( currentSeries->CY.cend()[-minSize] );
                for( auto helper : extra ) {
                    x1.push_back( helper->currentSeries->CX.cend()[-minSize] );
                    y1.push_back( helper->currentSeries->CY.cend()[-minSize] );
                }
                std::vector<double> x2;
                std::vector<double> y2;
                x2.push_back( currentSeries->CX.back() );
                y2.push_back( currentSeries->CY.back() );
                for( auto helper : extra ) {
                    x2.push_back( helper->currentSeries->CX.back() );
                    y2.push_back( helper->currentSeries->CY.back() );
                }

                CMatrix<double> Ax( 3, 1 );
                CMatrix<double> Ay( 3, 1 );
                if( LeastSquaresAffineTransform( Ax, Ay, x1, y1, x2, y2 ) ) {
                    double p[2];
                    if( SolveSystemOfTwoLinearEquations( p, Ax[0][0] - 1.0, Ax[1][0], Ax[2][0], Ay[0][0], Ay[1][0] - 1.0, Ay[2][0] ) ) {
                        PX = p[0];
                        PY = p[1];
                        qDebug() << "PX" << PX;
                        qDebug() << "PY" << PY;
                    }
                }
            }
        }
    }
    if( isGlobalPolarAlign && minSize > 1 ) {
        std::vector<StarXY> stars0;
        for( auto i = currentSeries->theDetectionResults.cend() - minSize; i < currentSeries->theDetectionResults.cend() - 1; i++ ) {
            const auto& dr = *i;
            if( dr.DetectionRegions.size() > 0 ) {
                stars0 = normalizedBrightStarsXY( dr, 50 );
                break;
            }
        }
        auto stars1 = normalizedBrightStarsXY( currentSeries->theDetectionResults.back(), 50 );
        Stars = align( stars0, stars1 );
        //currentSeries->DetectionResults.back().Image.reset();
        //currentSeries->DetectionResults.back().Mask.reset();

        std::vector<double> x1;
        std::vector<double> y1;
        std::vector<double> x2;
        std::vector<double> y2;
        for( const auto& s1 : stars1 ) {
            if( s1.MatchedIndex >= 0 ) {
                const auto& s0 = stars0[s1.MatchedIndex];
                x1.push_back( s0.X );
                y1.push_back( s0.Y );
                x2.push_back( s1.X );
                y2.push_back( s1.Y );
            }
        }
        if( x1.size() > 2 ) {
            CMatrix<double> Ax( 3, 1 );
            CMatrix<double> Ay( 3, 1 );
            if( LeastSquaresAffineTransform( Ax, Ay, x1, y1, x2, y2 ) ) {
                double p[2];
                if( SolveSystemOfTwoLinearEquations( p, Ax[0][0] - 1.0, Ax[1][0], Ax[2][0], Ay[0][0], Ay[1][0] - 1.0, Ay[2][0] ) ) {
                    PX1 = p[0];
                    PY1 = p[1];
                    qDebug() << "PX1" << PX1;
                    qDebug() << "PY1" << PY1;
                }
            }
        }
    }
}

DetectionResults CRawU16::DetectStars( int x0, int y0, int W, int H ) const
{
    auto image = GrayU16( x0, y0, W, H );

    const auto s = CRawU16::CalculateStatistics( image.get() ).stat( 0 );
    qDebug() << "Median:" << s.Median << "Sigma:" << s.Sigma;

    int th = s.Median + 6 * s.Sigma;
    int halfTh = s.Median + 3 * s.Sigma;


    DetectionResults results;
    auto& regions = results.DetectionRegions;
    auto mask = std::make_shared<CGrayImage>( W, H );
    for( int y = 0; y < H; y++ ) {
        const ushort* srcLine = image->ScanLine( y );
        uchar* dstLine = mask->ScanLine( y );
        for( int x = 0; x < W; x++ ) {
            const ushort* src = srcLine + x;
            uchar* dst = dstLine + x;

            if( dst[0] == 0 && src[0] >= th ) {
                auto region = std::make_shared<DetectionRegion>();
                int vmax, xmax, ymax, cmax, count;
                unsigned sum;
                std::tie( vmax, xmax, ymax, cmax, count, sum ) = starMaskWithMax( image.get(), x, y, halfTh, image->Width(), image->Height(), mask, 32 );
                region->Vmax = vmax;
                region->Xmax = xmax;
                region->Ymax = ymax;
                region->Cmax = cmax;
                region->AreaAtHalfDetectionThreshold = count;
                region->FluxAtHalfDetectionThreshold = sum;
                region->FluxAtHalfDetectionThreshold -= count * s.Median;
                region->Background = s.Median;
                int base = s.Median + ( vmax - s.Median ) / 10;
                if( base > halfTh ) {
                    starMask( image.get(), xmax, ymax, base, 64, mask );
                }
                starMask( image.get(), xmax, ymax, s.Median + ( vmax - s.Median ) / 2, 128, mask );
                starMask( image.get(), xmax, ymax, s.Median + 9 * ( vmax - s.Median ) / 10, 255, mask );
                regions.push_back( region );
            }
        }
    }
    results.Mask = mask;
    results.Image = image;
    return std::move( results );
}
