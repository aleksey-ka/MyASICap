// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include <locale>
#include <codecvt>

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "image.h"

static void fwprintf_no_trailing_zeroes( FILE* file, const wchar_t* name, float value )
{
    char buf[50];
    sprintf( buf, "%f", value );
    int pos = strlen( buf );
    for( int i = pos - 1; i > 0; i-- ) {
        if( buf[i] == '.' ) {
            buf[i] = 0;
            break;
        }
        if( buf[i] != '0' ) {
            buf[i + 1] = 0;
            break;
        }
    }
    fwprintf( file, L"%S %s\n", name, buf );
}

static void replace( std::string& str, const std::string& from, const std::string& to )
{
    size_t start_pos = str.find( from );
    if( start_pos != std::string::npos ) {
        str.replace( start_pos, from.length(), to );
    }
}

static std::string toString( std::wstring s )
{
    static std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converter;
    return converter.to_bytes( s );
}

std::shared_ptr<const Raw16Image> Raw16Image::LoadFromFile( const char* filePath )
{
    std::string infoFilePath( filePath );
    replace( infoFilePath, ".cfa", ".info" );

    std::map<std::wstring, std::wstring> map;

    std::wifstream info( infoFilePath );
    std::wstring line;
    while( std::getline( info, line ) ) {
        size_t pos = line.find_first_of( L" \t" );
        assert( pos != std::wstring::npos );
        map.insert( std::pair<std::wstring, std::wstring>(
            line.substr( 0, pos ),
            line.substr( pos + 1 ) ) );
    }
    info.close();

    ImageInfo imageInfo;

    int bits_per_pixel;
    swscanf( map[L"BITS_PER_PIXEL"].c_str(), L"%d", &bits_per_pixel );
    assert( bits_per_pixel == 16 );
    swscanf( map[L"WIDTH"].c_str(), L"%d", &imageInfo.Width );
    swscanf( map[L"HEIGHT"].c_str(), L"%d", &imageInfo.Height );
    swscanf( map[L"GAIN"].c_str(), L"%d", &imageInfo.Gain );
    swscanf( map[L"OFFSET"].c_str(), L"%d", &imageInfo.Offset );
    float exposure;
    swscanf( map[L"EXPOSURE"].c_str(), L"%f", &exposure );
    imageInfo.Exposure = (int)( 1000000 * exposure );
    swscanf( map[L"TEMPERATURE"].c_str(), L"%lf", &imageInfo.Temperature );
    imageInfo.Exposure = (int)( 1000000 * exposure );
    swscanf( map[L"TIMESTAMP"].c_str(), L"%I64d", &imageInfo.Timestamp );
    imageInfo.Camera = toString( map[L"CAMERA"] );

    auto result = std::make_shared<Raw16Image>( imageInfo );
    unsigned short* buf = (unsigned short*)result->Buffer();

    FILE* file = fopen( filePath, "rb" );
    fread( buf, 1, result->BufferSize(), file );
    fclose( file );

    for( int i = 0; i < result->BufferSize() / sizeof( buf[0] ); i++ ) {
       buf[i] >>= 4;
    }

    return result;
}

void Raw16Image::SaveToFile( const char* filePath ) const
{
    std::string infoFilePath( filePath );
    replace( infoFilePath, ".cfa", ".info" );
    FILE* info = fopen( infoFilePath.c_str(), "wt" );
    fwprintf( info, L"WIDTH %d\n", imageInfo.Width );
    fwprintf( info, L"HEIGHT %d\n", imageInfo.Height );
    fwprintf( info, L"BITS_PER_PIXEL 16\n" );
    fwprintf( info, L"OFFSET %d\n", imageInfo.Offset );
    fwprintf( info, L"PATTERN RGGB\n" );
    fwprintf( info, L"GAIN %d\n", imageInfo.Gain );
    fwprintf_no_trailing_zeroes( info, L"EXPOSURE", 0.000001f * imageInfo.Exposure );
    fwprintf_no_trailing_zeroes( info, L"TEMPERATURE", imageInfo.Temperature );
    fwprintf( info, L"TIMESTAMP %I64d\n", imageInfo.Timestamp );
    fwprintf( info, L"CAMERA %s\n", imageInfo.Camera.c_str() );
    fclose( info );

    FILE* out = fopen( filePath, "wb" );
    fwrite( Buffer(), 1, BufferSize(), out );
    fclose( out );
}
