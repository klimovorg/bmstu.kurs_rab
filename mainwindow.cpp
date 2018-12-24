#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include "QStandardItemModel"
#include "QStandardItem"
#include <cmath>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->textEdit_data_out->setText("Нажмите на кнопку 'Выполнить расчет'");
}

MainWindow::~MainWindow()
{
    delete ui;
}


//*******************************************************************************
//  фоновый поток
//*******************************************************************************
void MainWindow::on_pushButton_calc_clicked()
{
    ui->textEdit_data_out->setText("");

    int
            i = 1;
    double
            Tk = 0,
            Tcpu = 0,
            Tdisk = 0,
            Tsikle = 0,
            alpha_f = 0,
            delta_my = 0;
    double
            alpha_f1 = 0,
            betta = 1/(1- ui->doubleSpinBox_Y->value()),
            P = 1/ui->doubleSpinBox_M->value();

    alpha_f1 = ui->doubleSpinBox_K1->value()*
                minim(1/(2*ui->doubleSpinBox_t_k->value()),
                     ui->doubleSpinBox_C->value()/(betta*ui->doubleSpinBox_t_cpu->value()),
                     1/(betta * P * ui->doubleSpinBox_t_disk->value()))
                *((ui->doubleSpinBox_N->value()-1)/ui->doubleSpinBox_N->value());
    bool flag = true;
    while (flag)
    {
        Tk = (2*ui->doubleSpinBox_t_k->value()/(1-alpha_f1*ui->doubleSpinBox_t_k->value())),
        Tcpu = (betta*ui->doubleSpinBox_t_cpu->value())/
                (1-(betta*alpha_f1*pow(ui->doubleSpinBox_t_cpu->value()/ui->doubleSpinBox_C->value(),ui->doubleSpinBox_C->value()))),
        Tdisk = (betta*ui->doubleSpinBox_t_disk->value())/
                (1 - betta*P*alpha_f1*ui->doubleSpinBox_t_disk->value()),
        Tsikle =    ui->doubleSpinBox_T_dob->value()+
                    ui->doubleSpinBox_T_f->value()+
                    Tk+
                    Tcpu+
                    Tdisk,
        alpha_f = (ui->doubleSpinBox_N->value()-1)/Tsikle;

        delta_my = fabs(alpha_f1-alpha_f)/alpha_f;


        if (delta_my<ui->doubleSpinBox_delta->value())
        {
            flag = false;
        }else
        {
            double
                    S1 = fabs(alpha_f1-alpha_f)/ui->doubleSpinBox_K2->value();
            alpha_f1 = alpha_f1 - S1;
            ++i;
        }
    }

    ui->textEdit_data_out->append("----------Начало---------");
    ui->textEdit_data_out->append("Число итераций: " + QString::number(i));
    ui->textEdit_data_out->append("alpha_f1 = " + QString::number(alpha_f1));
    ui->textEdit_data_out->append("alpha_f = " + QString::number(alpha_f));
    ui->textEdit_data_out->append("Тр = " + QString::number(Tsikle - ui->doubleSpinBox_T_f->value()));
    ui->textEdit_data_out->append("Tk = " + QString::number(Tk));
    ui->textEdit_data_out->append("Tcpu = " + QString::number(Tcpu));
    ui->textEdit_data_out->append("Tdisk = " + QString::number(Tdisk));
    ui->textEdit_data_out->append("Tц = " + QString::number(Tsikle));
    double alpha = ui->doubleSpinBox_N->value()/Tsikle;
    ui->textEdit_data_out->append("Загрузка канала = " + QString::number(2*alpha*ui->doubleSpinBox_t_k->value()));
    ui->textEdit_data_out->append("Загрузка ЦП = " + QString::number((betta*alpha*ui->doubleSpinBox_t_cpu->value())/(ui->doubleSpinBox_C->value())));
    ui->textEdit_data_out->append("Загрузка дисков = " + QString::number(betta*P*alpha*ui->doubleSpinBox_t_disk->value()));
    ui->textEdit_data_out->append("Загрузка пользователя = " + QString::number(ui->doubleSpinBox_T_f->value()/Tsikle));
    ui->textEdit_data_out->append("Загрузка рабочей станции = " + QString::number((Tdisk+ui->doubleSpinBox_T_f->value())/Tsikle));
    ui->textEdit_data_out->append("----------Конец---------");
}

double MainWindow::minim(double x1,double x2,double x3)
{
    if (x1<x2)
    {
        if (x1<x3)
        {
            return x1;
        }else {
            return x3;
        }
    }else {
        if (x2<x3)
        {
            return x2;
        }else
        {
            return x3;
        }
    }
}


//*******************************************************************************
//  Модель ремонтника
//*******************************************************************************
void MainWindow::on_StartModel_clicked()
{
    double
            P0_1, P0_2, P0_3;
    long double
            mu_no,
            mu_o,
            psi,
            tno = ui->tno->text().toDouble();
    int
            c1,
            c2,
            c3,
            c=3,
            N;

    c1 = ui->C1->text().toInt();
    c2 = ui->C2->text().toInt();
    c3 = ui->C3->text().toInt();

    N = ui->N->text().toInt();

    mu_no = 1/ui->tno->text().toDouble();
    mu_o = 1/ui->to->text().toDouble();
    psi = mu_no/mu_o;

    //
    //      ---------------P0 start
    //

    P0_1 = findP(N, psi, 0, c1);
    P0_2 = findP(N, psi, 0, c2);
    P0_3 = findP(N, psi, 0, c3);
    //P0_3 = pow((a2 + b2), -1);

    //
    //      ----------------- Q
    //
    double Q_1 = 0;
    int C = c1;
    for (int k = C; k <= N; k++)
    {
        Q_1 = Q_1 + (k - C) * findP(N, psi, k, C);
    }

    double Q_2 = 0;
    C = c2;
    for (int k = C; k <= N; k++)
    {
        Q_2 = Q_2 + (k - C) * findP(N, psi, k, C);
    }

    double Q_3 = 0;
    C = c3;
    for (int k = C; k <= N; k++)
    {
        Q_3 = Q_3 + (k - C) * findP(N, psi, k, C);
    }

    //
    //      ----------------- L
    //
    double L_1 = 0;
    C = c1;
    for (int k = 1; k <= N; k++)
    {
        L_1 = L_1 + k * findP(N, psi, k, C);
    }

    double L_2 = 0;
    C = c2;
    for (int k = 1; k <= N; k++)
    {
        L_2 = L_2 + k * findP(N, psi, k, C);
    }

    double L_3 = 0;
    C = c3;
    for (int k = 1; k <= N; k++)
    {
        L_3 = L_3 + k * findP(N, psi, k, C);
    }

    //расчет U
    double U_1 = L_1 - Q_1;
    double U_2 = L_2 - Q_2;
    double U_3 = L_3 - Q_3;
    //расчет r0
    C = c1;
    double r0_1 = U_1 / C;
    C = c2;
    double r0_2 = U_2 / C;
    C = c3;
    double r0_3 = U_3 / C;
    //расчет Tp
    double Tp_1 = L_1 * tno / (N - L_1);
    double Tp_2 = L_2 * tno / (N - L_2);
    double Tp_3 = L_3 * tno / (N - L_3);
    //расчет W
    double to = ui->to->text().toDouble();
    double W_1 = Tp_1 - to;
    double W_2 = Tp_2 - to;
    double W_3 = Tp_3 - to;
    //расчет Tcircle
    double Tcircle_1 = Tp_1 + tno;
    double Tcircle_2 = Tp_2 + tno;
    double Tcircle_3 = Tp_3 + tno;
    //расчет n
    double n_1 = N - L_1;
    double n_2 = N - L_2;
    double n_3 = N - L_3;
    //расчет re
//    double re_1 = to / Tcircle_1;
//    double re_2 = to / Tcircle_2;
//    double re_3 = to / Tcircle_3;
    double re_1 = n_1 / N;
    double re_2 = n_2 / N;
    double re_3 = n_3 / N;
    //расчет re/r0
    double rero_1 = re_1 / r0_1;
    double rero_2 = re_2 / r0_2;
    double rero_3 = re_3 / r0_3;
    //расчет Y
    double S1 = ui->S1->text().toDouble(),
            S = ui->S_rem->text().toDouble(),
            Y_1 = (c1) * S1 + L_1 * S;
    double Y_2 = (c2) * S1 + L_2 * S;
    double Y_3 = (c3) * S1 + L_3 * S;

    //Выгрузка
    QStandardItemModel *model = new QStandardItemModel;
    QStandardItem *item;

    //Заголовки столбцов
    QStringList horizontalHeader;
    horizontalHeader.append("Вар 1");
    horizontalHeader.append("Вар 2");
    horizontalHeader.append("Вар 3");
    horizontalHeader.append("Информация");

    //Заголовки строк
    QStringList verticalHeader;
    verticalHeader.append("с");
    verticalHeader.append("P0");
    verticalHeader.append("Q");
    verticalHeader.append("L");
    verticalHeader.append("U");
    verticalHeader.append("p0");
    verticalHeader.append("n");
    verticalHeader.append("pe");
    verticalHeader.append("W");
    verticalHeader.append("Tp");
    verticalHeader.append("Tц");
    verticalHeader.append("pe/p0");
    verticalHeader.append("Y");

    model->setHorizontalHeaderLabels(horizontalHeader);
    model->setVerticalHeaderLabels(verticalHeader);

    //c
    item = new QStandardItem(QString::number(c1));
    model->setItem(0, 0, item);
    item = new QStandardItem(QString::number(c2));
    model->setItem(0, 1, item);
    item = new QStandardItem(QString::number(c3));
    model->setItem(0, 2, item);
    item = new QStandardItem(QString("кол-во ремонтников"));
    model->setItem(0, 3, item);

    //P0
    item = new QStandardItem(QString::number(P0_1));
    model->setItem(1, 0, item);
    item = new QStandardItem(QString::number(P0_2));
    model->setItem(1, 1, item);
    item = new QStandardItem(QString::number(P0_3));
    model->setItem(1, 2, item);
    item = new QStandardItem(QString("Вероятность состояния 0"));
    model->setItem(1, 3, item);

    //Q
    item = new QStandardItem(QString::number(Q_1));
    model->setItem(2, 0, item);
    item = new QStandardItem(QString::number(Q_2));
    model->setItem(2, 1, item);
    item = new QStandardItem(QString::number(Q_3));
    model->setItem(2, 2, item);
    item = new QStandardItem(QString("среднее кол-во ПК в очереди на ремонт"));
    model->setItem(2, 3, item);

    //L
    item = new QStandardItem(QString::number(L_1));
    model->setItem(3, 0, item);
    item = new QStandardItem(QString::number(L_2));
    model->setItem(3, 1, item);
    item = new QStandardItem(QString::number(L_3));
    model->setItem(3, 2, item);
    item = new QStandardItem(QString("среднее кол-во ПК в неисправном состонии"));
    model->setItem(3, 3, item);

    //U
    item = new QStandardItem(QString::number(U_1));
    model->setItem(4, 0, item);
    item = new QStandardItem(QString::number(U_2));
    model->setItem(4, 1, item);
    item = new QStandardItem(QString::number(U_3));
    model->setItem(4, 2, item);
    item = new QStandardItem(QString("среднее кол-во ПК которое непосредственно ремонтируется специалистами"));
    model->setItem(4, 3, item);

    //p0
    item = new QStandardItem(QString::number(r0_1));
    model->setItem(5, 0, item);
    item = new QStandardItem(QString::number(r0_2));
    model->setItem(5, 1, item);
    item = new QStandardItem(QString::number(r0_3));
    model->setItem(5, 2, item);
    item = new QStandardItem(QString("коэф загрузки одного специалиста"));
    model->setItem(5, 3, item);

    //n
    item = new QStandardItem(QString::number(n_1));
    model->setItem(6, 0, item);
    item = new QStandardItem(QString::number(n_2));
    model->setItem(6, 1, item);
    item = new QStandardItem(QString::number(n_3));
    model->setItem(6, 2, item);
    item = new QStandardItem(QString("среднее кол-во исправных компьютеров"));
    model->setItem(6, 3, item);

    //pe
    item = new QStandardItem(QString::number(re_1));
    model->setItem(7, 0, item);
    item = new QStandardItem(QString::number(re_2));
    model->setItem(7, 1, item);
    item = new QStandardItem(QString::number(re_3));
    model->setItem(7, 2, item);
    item = new QStandardItem(QString("коэф загрузки ПК"));
    model->setItem(7, 3, item);

    //W
    item = new QStandardItem(QString::number(W_1));
    model->setItem(8, 0, item);
    item = new QStandardItem(QString::number(W_2));
    model->setItem(8, 1, item);
    item = new QStandardItem(QString::number(W_3));
    model->setItem(8, 2, item);
    item = new QStandardItem(QString("среднее время нахождения ПК на ремонт"));
    model->setItem(8, 3, item);

    //Tp
    item = new QStandardItem(QString::number(Tp_1));
    model->setItem(9, 0, item);
    item = new QStandardItem(QString::number(Tp_2));
    model->setItem(9, 1, item);
    item = new QStandardItem(QString::number(Tp_3));
    model->setItem(9, 2, item);
    item = new QStandardItem(QString("среднее время пребывания ПК в неисправном состоянии (очередь + ремонт)"));
    model->setItem(9, 3, item);

    //Tц
    item = new QStandardItem(QString::number(Tcircle_1));
    model->setItem(10, 0, item);
    item = new QStandardItem(QString::number(Tcircle_2));
    model->setItem(10, 1, item);
    item = new QStandardItem(QString::number(Tcircle_3));
    model->setItem(10, 2, item);
    item = new QStandardItem(QString("среднее время цикла для ПК"));
    model->setItem(10, 3, item);

    //pe/p0
    item = new QStandardItem(QString::number(rero_1));
    model->setItem(11, 0, item);
    item = new QStandardItem(QString::number(rero_2));
    model->setItem(11, 1, item);
    item = new QStandardItem(QString::number(rero_3));
    model->setItem(11, 2, item);
    item = new QStandardItem(QString("если = 1, то система сбалансирована"));
    model->setItem(11, 3, item);

    //Y
    item = new QStandardItem(QString::number(Y_1));
    model->setItem(12, 0, item);
    item = new QStandardItem(QString::number(Y_2));
    model->setItem(12, 1, item);
    item = new QStandardItem(QString::number(Y_3));
    model->setItem(12, 2, item);
    item = new QStandardItem(QString("Убытки организации"));
    model->setItem(12, 3, item);




    ui->tableView_out->setModel(model);

    ui->tableView_out->resizeRowsToContents();
    ui->tableView_out->resizeColumnsToContents();


}


long double MainWindow::factorial(int N)
{
    if(N < 0) // если пользователь ввел отрицательное число
        return 0; // возвращаем ноль
    if (N == 0) // если пользователь ввел ноль,
        return 1; // возвращаем факториал от нуля - не удивляетесь, но это 1 =)
    else // Во всех остальных случаях
        return N * factorial(N - 1); // делаем рекурсию.
}

double MainWindow::findP(int N,double psi,int k,int c)
{
                double a = 0, b = 0;
                double a1 = 0, b1 = 0, a2 = 0, b2 = 0;
                //расчет P0. Если будет сильно тупить - перекинуть в отдельную функцию и сохранать значение а не пересчитывать
                int k1 = 0;
                while (k1 <= c)
                {
                    a = (factorial(N) / (factorial(k1) * factorial(N - k1)));
                    a1 = (double)a;
                    a1 =  a1 * pow(psi, k1);
                    a2 = a2 + a1;
                    k1++;
                }

                while (k1 <= N)
                {
                    b = factorial(N)  / ( factorial(c) * factorial(N - k1));
                    b1 = (double)b;
                    b1 = b1 * pow(psi, k1) / pow(c, (k1 - c));
                    b2 = b2 + b1;
                    k1++;
                }

                double P0 = pow((a2 + b2), -1);
                // собственно расчет вероятности состояния
                if (k == 0)
                {
                    return P0;
                }
                if (k <= c)
                {
                    double pkc = factorial(N) / (factorial(k) * factorial(N - k));
                    double pkc1 = (double)pkc * pow(psi, k) * P0;
                    return pkc1;
                }
                if (k>c)
                {
                    double pck = factorial(N) / (factorial(c) * factorial(N - k));
                    double pck1 = (double)pck * pow(psi, k) * P0 / pow(c, (k - c));
                    return pck1;
                }
                return -1;
}
