#ifndef FOCUSERCONTROLS_H
#define FOCUSERCONTROLS_H

#include <QFrame>

namespace Ui {
class FocuserControls;
}

class FocuserControls : public QFrame
{
    Q_OBJECT

public:
    explicit FocuserControls(QWidget *parent = nullptr);
    ~FocuserControls();

private:
    Ui::FocuserControls *ui;
};

#endif // FOCUSERCONTROLS_H
