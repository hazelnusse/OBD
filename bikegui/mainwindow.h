#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class QAction;

class MainWindow : public QMainWindow
{
  Q_OBJECT;

public:
  MainWindow();

private slots:
  void about(void);

private:
  // Member functions
  void createActions(void);
  void createMenus(void);

  // Menu's
  QMenu *fileMenu;
  QMenu *helpMenu;

  // Action's
  QAction *newAction;
  QAction *openAction;
  QAction *saveAction;
  QAction *saveAsAction;
  QAction *quitAction;
  QAction *aboutAction;
  QAction *aboutQtAction;
};

#endif
