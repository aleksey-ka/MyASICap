// Copyright (C) 2021 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#ifndef IMAGE_DEBAYER_H
#define IMAGE_DEBAYER_H

#include <limits.h>
#include <vector>

class CDebayer_RawU16 {
public:
    CDebayer_RawU16( const unsigned short* _raw, int _width, int _height ) :
        raw( _raw ), width( _width ), height( _height )
    {
    }

    unsigned short MaxValue = 0;
    unsigned MaxCount = 0;
    unsigned short MinValue = USHRT_MAX;
    unsigned MinCount = 0;

protected:
    const unsigned short* raw;
    int width;
    int height;

    // Fast statistics (calculated for each pixel on each frame)
    unsigned short addToStatistics( unsigned short value );
};

inline unsigned short CDebayer_RawU16::addToStatistics( unsigned short value )
{
    if( value >= MaxValue ) {
        if( value == MaxValue ) {
            MaxCount++;
        } else {
            MaxValue = value;
            MaxCount = 1;
        }
    } else if( value <= MinValue ) {
        if( value == MinValue ) {
            MinCount++;
        } else {
            MinValue = value;
            MinCount = 1;
        }
    }
    return value;
}

#define CFA_CHANNEL_AT( x, y ) ( x % 2 | y % 2 << 1 )

#endif // IMAGE_DEBAYER_H
