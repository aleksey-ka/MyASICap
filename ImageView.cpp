// Copyright (C) 2020 Aleksey Kalyuzhny. Released under the terms of the
// GNU General Public License version 3. See <http://www.gnu.org/licenses/>

#include "ImageView.h"

ImageView::ImageView( QWidget *parent ) : QLabel( parent )
{
}

void ImageView::mousePressEvent( QMouseEvent* event )
{
    if( pixmap() != 0 ) {
        auto cr = contentsRect();
        auto m = margin();
        cr.adjust( m,  m, -m, -m );
        auto rect = style()->itemPixmapRect( cr, alignment(), *pixmap() );
        int x = event->x() - rect.x();
        int y = event->y() - rect.y();
        if( x >= 0 && y >= 0 && x < rect.width() && y < rect.height() ) {
            emit imagePressed( x, y, event->button(), event->modifiers() );
        }
    }
}
