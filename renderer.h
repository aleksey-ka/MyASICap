#ifndef RENDERER_H
#define RENDERER_H

#include <QImage>
#include <QPixmap>

class Renderer {
public:
    Renderer( const ushort* raw, int width, int height );

    QPixmap RenderHalfResWithHistogram();
    QPixmap RenderHistogram();

    QPixmap RenderRectHalfRes( int cx, int cy, int W, int H );

private:
    const ushort* raw;
    int width;
    int height;

    std::vector<uint> histR;
    std::vector<uint> histG;
    std::vector<uint> histB;
};

#endif // RENDERER_H
