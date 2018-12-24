#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_calc_clicked();

    void on_StartModel_clicked();

private:
    double minim(double x1,double x2,double x3);
    Ui::MainWindow *ui;
    //факториал числа
    long double factorial(int N);
    double findP(int N,double psi,int k,int c);
};

#endif // MAINWINDOW_H
