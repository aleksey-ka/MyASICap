// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#ifndef RENDERER_H
#define RENDERER_H

#include <QImage>
#include <QPixmap>

enum TRenderingMethod {
    RM_HalfResolution,
    RM_FullResolution
};

class Renderer {
public:
    Renderer( const ushort* raw, int width, int height );

    QPixmap Render( TRenderingMethod, int x = 0, int y = 0, int w = 0, int h = 0 );
    QPixmap RenderHistogram();

    QPixmap RenderGrayScale();
private:
    const ushort* raw;
    int width;
    int height;

    // Histogram data
    std::vector<uint> histR;
    std::vector<uint> histG;
    std::vector<uint> histB;

    ushort maxValue = 0;
    uint maxCount = 0;
    ushort minValue = USHRT_MAX;
    uint minCount = 0;
};

#endif // RENDERER_H
