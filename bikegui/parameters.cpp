#include <QtGui>
#include "parameters.h"

WhippleParameter::WhippleParameter(QWidget *parent)
  : QWidget(parent)
{
  comboBox = new QComboBox(this);
  comboBox->addItem(tr("Gyrostat parameters"));
  comboBox->addItem(tr("Franke parameters"));
  comboBox->addItem(tr("Benchmark parameters"));

  paramBox = new QGroupBox(this);

  layout = new QVBoxLayout(this);
  layout->addWidget(comboBox);
  layout->addWidget(paramBox);
}
