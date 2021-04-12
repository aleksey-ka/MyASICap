// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#ifndef IMAGE_RGBIMAGE_H
#define IMAGE_RGBIMAGE_H

#include <vector>

template<typename T, int numOfChannels = 1>
class CPixelBuffer {
public:
    CPixelBuffer( int width, int height );

    int Width() const { return width; }
    int Height() const { return height; }

    const T* Pixels() const { return buffer.data(); };
    T* Pixels() { return buffer.data(); };

    const T* ScanLine( int y ) const { return buffer.data() + y * stride; }
    T* ScanLine( int y ) { return buffer.data() + y * stride; }

    const T* Ptr( int x, int y ) const { return ScanLine( y ) + numOfChannels * x; }
    T* Ptr( int x, int y ) { return ScanLine( y ) + numOfChannels * x; }

    const T& At( int x, int y, int ch = 0 ) const { return ( ScanLine( y ) + numOfChannels * x )[ch]; }
    T& At( int x, int y, int ch = 0 ) { return ( ScanLine( y ) + numOfChannels * x )[ch]; }

    size_t Count() { return buffer.size(); }

protected:
    int width;
    int height;
    int stride;
    std::vector<T> buffer;
};

template<typename T, int numOfChannels>
inline CPixelBuffer<T, numOfChannels>::CPixelBuffer( int _width, int _height ) :
    width( _width ), height( _height )
{
    stride = numOfChannels * width;
    buffer.resize( sizeof( T ) * stride * height );
}

class CRgbU16Image : public CPixelBuffer<unsigned short, 3> {
public:
    using CPixelBuffer::CPixelBuffer;

    const unsigned short* RgbPixels() const { return buffer.data(); };
    unsigned short* RgbPixels() { return buffer.data(); };

    int Stride() const { return stride; }
};

class CRgbImage : public CPixelBuffer<unsigned char, 3> {
public:
    using CPixelBuffer::CPixelBuffer;

    const unsigned char* RgbPixels() const { return buffer.data(); };
    unsigned char* RgbPixels() { return buffer.data(); };

    int ByteWidth() const { return stride; }
};

class CGrayU16Image : public CPixelBuffer<unsigned short, 1> {
public:
    using CPixelBuffer::CPixelBuffer;

    int Stride() const { return stride; }
};

class CGrayImage : public CPixelBuffer<unsigned char, 1> {
public:
    using CPixelBuffer::CPixelBuffer;

    int ByteWidth() const { return stride; }
};

#endif // IMAGE_RGBIMAGE_H
