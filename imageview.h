#ifndef IMAGEVIEW_H
#define IMAGEVIEW_H

#include <QLabel>
#include <QMouseEvent>
#include <QStyle>
#include "image.h"

class ImageView : public QLabel
{
    Q_OBJECT
public:
    explicit ImageView(QWidget *parent = nullptr);

    void SetImage( std::shared_ptr<const Raw16Image> _image ) { image = _image; }

protected:
    void mousePressEvent( QMouseEvent *event ) override
    {
        if( event->button() == Qt::LeftButton ) {
            auto modifiers = event->modifiers();
            if( modifiers.testFlag( Qt::ControlModifier ) && modifiers.testFlag( Qt::ShiftModifier ) ) {
                auto cr = contentsRect();
                auto m = margin();
                cr.adjust( m,  m, -m, -m );
                auto rect = style()->itemPixmapRect( cr, alignment(), *pixmap() );
                render( event->x() - rect.x(), event->y() - rect.y() );
            }
        }
    }

signals:

private:
    std::shared_ptr<const Raw16Image> image;
    void render( int x, int y );
};

#endif // IMAGEVIEW_H
