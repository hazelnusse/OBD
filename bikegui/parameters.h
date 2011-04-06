#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <QWidget>

class QVBoxLayout;
class QComboBox;
class QGroupBox;

class WhippleParameter : public QWidget
{
  Q_OBJECT

public:
  WhippleParameter(QWidget *parent = 0);

private:
  QVBoxLayout *layout;
  QComboBox *comboBox;
  QGroupBox *paramBox;
};
#endif
