#include "focusercontrols.h"
#include "ui_focusercontrols.h"

FocuserControls::FocuserControls(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::FocuserControls)
{
    ui->setupUi(this);
}

FocuserControls::~FocuserControls()
{
    delete ui;
}
