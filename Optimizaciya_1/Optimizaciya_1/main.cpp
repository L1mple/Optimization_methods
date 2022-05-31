#include <iostream>                              
#include <cmath>   



using namespace std;

extern double Name(4.0);
extern double Surname(5.0);
extern double Patronymic(10.0);
extern double Count(1);


double F1(double x) {
    return Surname * x;
}

double F2(double x) {
    return Name * x * (-1);
}

double F3(double x) {
    return ((Name) * (x * x) + Patronymic * x + Surname - x * x * x);
}
double MainFunction(double x) {
    Count += 1;
    return exp(F1(x)) + exp(F2(x)) + F3(x);
}
int main() {
    double epsilon = 0.001;
    double minimum = 100;
    double a = -100;           //левая граница интервала 
    double b = 100;            //правая граница интервала
    double delta = 0.0001;
    double stopper = (b - a) / 2; // Критерий остановки, как только он меньше заданной погрешности метод дихотомии останавливается
    double c, d;
    double passive;
    for (double i = -100; i <= 100; i = i + 0.001) {
        if (MainFunction(i) < minimum) {
            minimum = MainFunction(i);
            passive = i;
        }

    };
    cout << "Min of passive seacrh = " << passive << endl;

    while (stopper > epsilon) {
        c = (a + b) / 2 - delta / 2;
        d = (a + b) / 2 + delta / 2;
        if (MainFunction(c) <= MainFunction(d)) {
            a = a;
            b = d;
        };
        if (MainFunction(c) > MainFunction(d)) {
            a = c;
            b = b;
        };
        stopper = (b - a) / 2;
    };
    cout << "Min of dichotomy = " << (a + b) / 2 << endl;
    cout << MainFunction(passive) << endl;
    cout << MainFunction((a + b) / 2) << endl;
    cout << "Number of function calls = " << Count << endl;


    return 0;
}