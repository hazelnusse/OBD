#include <QtGui>
#include "mainwindow.h"
#include "OBDConfig.h"
#include <iostream>

MainWindow::MainWindow()
{
  createActions();
  createMenus();
}

void MainWindow::about(void)
{
  QString verstring = QString(tr(
        "<h2>Open Bicycle Dynamics</h2>"
        "Version %1.%2 commit %3"
        "<p>Copyright &copy; 2010-2011 Dale Lukas Peterson"
        "<p><a href=\"mailto:%4\">%4</a>"))
                    .arg(OBD_VERSION_MAJOR)
                    .arg(OBD_VERSION_MINOR)
                    .arg(OBD_VERSION_COMMIT)
                    .arg(OBD_AUTHOR_CONTACT);
  QMessageBox::about(this, tr("About OBD"), verstring);
}

void MainWindow::createActions(void)
{
  // New action
  newAction = new QAction(tr("&New"), this);
  newAction->setShortcut(QKeySequence::New);
  newAction->setStatusTip(tr("Create a new set of bicycle parameters"));
  // connect(newAction, SIGNAL(triggered()), this, SLOT(newFile()));

  // Open action
  openAction = new QAction(tr("&Open"), this);
  openAction->setShortcut(QKeySequence::Open);
  openAction->setStatusTip(tr("Open a set of bicycle parameters"));
  // connect(openAction, SIGNAL(triggered()), this, SLOT(about()));

  // Save action
  saveAction = new QAction(tr("&Save"), this);
  saveAction->setShortcut(QKeySequence::Save);
  saveAction->setStatusTip(tr("Save bicycle parameters to file"));
  // connect(saveAction, SIGNAL(triggered()), this, SLOT(save()));

  // Save as action
  saveAsAction = new QAction(tr("Save &As"), this);
  saveAsAction->setShortcut(QKeySequence::SaveAs);
  saveAsAction->setStatusTip(tr("Save bicycle parameters to file "));
  // connect(saveAsAction, SIGNAL(triggered()), this, SLOT(saveAs()));

  // Quit action
  quitAction = new QAction(tr("&Quit"), this);
  quitAction->setShortcut(QKeySequence::Quit);
  quitAction->setStatusTip(tr("Quit"));
  connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));

  // About action
  aboutAction = new QAction(tr("&About"), this);
  aboutAction->setStatusTip(tr("About OBD"));
  connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));

  // About Qt action
  aboutQtAction = new QAction(tr("About &Qt"), this);
  aboutQtAction->setStatusTip(tr("About Qt"));
  connect(aboutQtAction, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void MainWindow::createMenus(void)
{
  // Add the file menu
  fileMenu = menuBar()->addMenu(tr("&File"));
  // Add actions to file menu
  fileMenu->addAction(newAction);
  fileMenu->addAction(openAction);
  fileMenu->addAction(saveAction);
  fileMenu->addAction(saveAsAction);
  fileMenu->addAction(quitAction);

  // Add the help menu
  helpMenu = menuBar()->addMenu(tr("&Help"));
  // Add actions to help menu
  helpMenu->addAction(aboutAction);
  helpMenu->addAction(aboutQtAction);
}
