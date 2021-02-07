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

class Pixels16BitUncompressed : public ImageFileFormat {
public:
    virtual std::shared_ptr<const Raw16Image> Load( const char* filePath, const ImageInfo& ) const override;
    virtual void Save( const char* filePath, const Raw16Image* ) const override;
};

std::shared_ptr<const Raw16Image> Pixels16BitUncompressed::Load( const char* filePath, const ImageInfo& imageInfo ) const
{
    auto result = std::make_shared<Raw16Image>( imageInfo );

    FILE* file = fopen( filePath, "rb" );
    fread( result->Buffer(), 1, result->BufferSize(), file );
    fclose( file );

    return result;
}

void Pixels16BitUncompressed::Save( const char* filePath, const Raw16Image* image ) const
{
    FILE* out = fopen( filePath, "wb" );
    fwrite( image->Buffer(), 1, image->BufferSize(), out );
    fclose( out );
}

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
    replace( infoFilePath, ".pixels", ".info" );

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

    assert( map[L"PIXEL_FORMAT"] == L"u16" );
    imageInfo.CFA = toString( map[L"PIXEL_CFA"] );
    swscanf( map[L"IMAGE_WIDTH"].c_str(), L"%d", &imageInfo.Width );
    swscanf( map[L"IMAGE_HEIGHT"].c_str(), L"%d", &imageInfo.Height );
    swscanf( map[L"CAMERA_GAIN"].c_str(), L"%d", &imageInfo.Gain );
    swscanf( map[L"CAMERA_OFFSET"].c_str(), L"%d", &imageInfo.Offset );
    float exposure;
    swscanf( map[L"CAMERA_EXPOSURE"].c_str(), L"%f", &exposure );
    imageInfo.Exposure = (int)( 1000000 * exposure );
    swscanf( map[L"CAMERA_TEMPERATURE"].c_str(), L"%lf", &imageInfo.Temperature );
    imageInfo.Exposure = (int)( 1000000 * exposure );
    swscanf( map[L"TIMESTAMP"].c_str(), L"%I64d", &imageInfo.Timestamp );
    swscanf( map[L"SERIES_ID"].c_str(), L"%I64d", &imageInfo.SeriesId );
    imageInfo.Camera = toString( map[L"CAMERA"] );
    imageInfo.Channel = toString( map[L"CHANNEL"] );
    imageInfo.FilterDescription = toString( map[L"FILTER"] );

    Pixels16BitUncompressed uncompressed;
    return uncompressed.Load( filePath, imageInfo );
}

void Raw16Image::SaveToFile( const char* filePath, const ImageFileFormat* fileFormat ) const
{
    std::string infoFilePath( filePath );
    auto pos = infoFilePath.rfind( '.' );
    infoFilePath.replace( pos, infoFilePath.length() - pos, ".info" );
    FILE* info = fopen( infoFilePath.c_str(), "wt" );
    fwprintf( info, L"IMAGE_WIDTH %d\n", imageInfo.Width );
    fwprintf( info, L"IMAGE_HEIGHT %d\n", imageInfo.Height );
    fwprintf( info, L"PIXEL_FORMAT u16\n" );
    if( not imageInfo.CFA.empty() ) {
        fwprintf( info, L"PIXEL_CFA %s\n", imageInfo.CFA.c_str() );
    }
    if( not imageInfo.Channel.empty() ) {
        fwprintf( info, L"CHANNEL %s\n", imageInfo.Channel.c_str() );
    }
    if( not imageInfo.FilterDescription.empty() ) {
        fwprintf( info, L"FILTER %s\n", imageInfo.FilterDescription.c_str() );
    }
    fwprintf( info, L"CAMERA %s\n", imageInfo.Camera.c_str() );
    fwprintf( info, L"CAMERA_GAIN %d\n", imageInfo.Gain );
    fwprintf( info, L"CAMERA_OFFSET %d\n", imageInfo.Offset );
    fwprintf_no_trailing_zeroes( info, L"CAMERA_EXPOSURE", 0.000001f * imageInfo.Exposure );
    fwprintf_no_trailing_zeroes( info, L"CAMERA_TEMPERATURE", imageInfo.Temperature );
    fwprintf( info, L"SERIES_ID %I64d\n", imageInfo.SeriesId );
    fwprintf( info, L"TIMESTAMP %I64d\n", imageInfo.Timestamp );
    fclose( info );

    if( fileFormat == 0 ) {
        static const Pixels16BitUncompressed uncompressed;
        fileFormat = &uncompressed;
    }

    fileFormat->Save( filePath, this );

    // Check saved file
    /*auto saved = fileFormat->Load( filePath, imageInfo );
    auto size = saved->BufferSize();
    auto src = Buffer();
    auto dst = saved->Buffer();
    assert( BufferSize() == size );
    for( int i = 0; i < size; i++ ) {
        assert( dst[i] == src[i] );
    }*/
}
