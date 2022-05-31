#include <iostream>                              
#include <cmath>   
#include <iomanip>


using namespace std;

static double Name(4.0);
static double Surname(5.0);
static double Patronymic(10.0);
static double Count(0);


double F1(double x) {
    return Surname * x;
}

double F2(double x) {
    return Name * x * (-1);
}

double F3(double x) {
    return ((Name) * (x*x) + Patronymic * x + Surname - x*x*x);
}
double MainFunction(double x) {
    Count += 1;
    return exp(F1(x)) + exp(F2(x)) + F3(x);
}
double MainFunctionDif(double x) {
    Count += 1;
    return 5 * exp(F1(x)) - 4 * exp(F2(x)) + 8*x + 10 - 3*x*x;
}
double MainFunctionDif2(double x) {
    Count += 1;
    return 25 * exp(F1(x)) + 16 * exp(F2(x)) + 8  - 6 * x;
}

// Fibonacci
int f(int n)
{
    int a = 0;
    int b = 1;
    for (int i = 0; i < n; i++)
    {
        a = a + b;
        b = a - b;
    }
    return a;
}


int main(){
    double epsilon = 0.000001;    
    double minimum = 100;
    double a = -100;           //
    double b = 100;            //
    double delta; 
    double stopper = (b - a) / 2; // Когда эта величина меньше заранее заданной погрешности, метод дихотоммии останавливается
    double c, d;
    double passive;
    double F, FC, FD;
    int iteration = 1;
    for (double i = -100; i <= 100; i = i + 0.001) {
        F = MainFunction(i);
        Count -= 1;
        if (F < minimum) {
            minimum = F;
            passive = i;
            
        }

    };
    
    cout << "Min of passive seacrh = " << passive << endl;
    cout << "The dichotomy method" << endl;
    cout << setw(9) << "Iteration " << setw(12) << "a "
        << setw(9) << "b" << setw(12) << "Stopper"
        << setw(9) << "c" << setw(12) << "F(c)"
        << setw(9) << "d" << setw(12) << "F(d)" << endl;
    while (stopper > epsilon) {
        delta = (b - a) / 4;
        c = (a + b) / 2 - delta / 2;
        d = (a + b) / 2 + delta / 2;
        FC = MainFunction(c);
        FD = MainFunction(d);
        cout << setw(10) << iteration << setw(13) << a
            << setw(10) << b << setw(13) << stopper
            << setw(10) << c << setw(13) << FC
            << setw(10) << d << setw(13) << FD << endl;
        if (FC <= FD) {
            a = a;
            b = d;
        };
        if (FC > FD) {
            a = c;
            b = b;
        };
        
        stopper = (b - a) / 2;
        iteration += 1;
    };
    cout << "Min of dichotomy = " << setprecision(4) << (a + b) / 2 << endl;
    //cout << MainFunction(passive) << endl;
    //cout << MainFunction((a + b) / 2) << endl;
    cout << "Number of function calls = " << Count << endl;
    Count = 0;
    cout << "Golden ratio method" << endl;
    a = -100;
    b = 100;
    iteration = 1;
    stopper = (b - a) / 2;
    FC = FD = 0;
    cout << setw(9) << "Iteration " << setw(12) << "a "
        << setw(9) << "b" << setw(12) << "Stopper"
        << setw(9) << "c" << setw(12) << "F(c)"
        << setw(9) << "d" << setw(12) << "F(d)" << endl;
    while (stopper > epsilon) {
        
        c = ((3 - sqrt(5)) / 2) * (b - a) + a;
        d = ((sqrt(5) - 1) / 2) * (b - a) + a;
        if (FC == 0) {
            FC = MainFunction(c);
        }
        if (FD == 0) {
            FD = MainFunction(d);
        }
        cout << setw(10) << setprecision(6) << iteration << setw(13) << setprecision(6) << a
            << setw(10) << setprecision(6) << b << setw(13) << setprecision(6) << stopper
            << setw(10) << setprecision(6) << c << setw(13) << setprecision(6) << FC
            << setw(10) << setprecision(6) << d << setw(13) << setprecision(6) << FD << endl;


        if (FC <= FD) {
            a = a;
            b = d;
            FD = FC;
            FC = 0;
        };
        if (FC > FD) {
            a = c;
            b = b;
            FC = FD;
            FD = 0;
        };


        stopper = (b - a) / 2;
        iteration += 1;
    };
    cout << "Min of golden ratio = " << setprecision(4) << (a + b) / 2 << endl;
    cout << "Number of function calls = " << Count << endl;
    Count = 0;
    cout << "Fibonacci method" << endl;
    a = -100;
    b = 100;
    iteration = 1;
    FC = FD = 0;
    int n = 1;
    
    while ((f(n+2)) < ((b-a)/epsilon)) {
        n += 1;
    }
    cout << "Suitable n for method: " << n << endl;
    //cout << "Suitable n for method: " << f(n + 1 - 1) << endl;
    //cout << "Suitable n for method: " << f(n + 3 - 1) << endl;
    cout << setw(12) << "Iteration " << setw(16) << "a "
        << setw(12) << "b" 
        << setw(12) << "c" << setw(16) << "F(c)"
        << setw(12) << "d" << setw(16) << "F(d)" << endl;
    
    for (int i = 1; i < (n+1); i++) {
        long double e, r, t;
        e = f(n + 1 - i);
        r = f(n + 2 - i);
        t = f(n + 3 - i);
        c = a + (b - a) * (e / t);
     
        d = a + (b - a) * (r / t);
        
        if (FC == 0) {
            FC = MainFunction(c);
            
        }
        if (FD == 0) {
            FD = MainFunction(d);
            
        }



        cout << setw(12) << setprecision(7) << iteration << setw(16) << setprecision(7) << a
            << setw(12) << setprecision(7) << b 
            << setw(12) << setprecision(7) << c << setw(16) << setprecision(7) << FC
            << setw(12) << setprecision(7) << d << setw(16) << setprecision(7) << FD << endl;

        if (FC <= FD) {
            a = a;
            b = d;
            FD = FC;
            FC = 0;
        };
        if (FC > FD) {
            a = c;
            b = b;
            FC = FD;
            FD = 0;
        };


        iteration += 1;
    }
    cout << "Min of Fibonacci = " << setprecision(4) << (a + b) / 2 << endl;
    cout << "Number of function calls = " << Count << endl;
    Count = 0;


    cout << "Tangent method" << endl;
    
    a = -100;
    b = 100;
    iteration = 1;
    FC = FD = 0;
    double FA, FB, DFA, DFB, DFC;
    FA = FB = DFA = DFB = DFC = 0;
    epsilon = 0.000001;;
    cout << setw(12) << "Iteration " << setw(16) << "a "
        << setw(12) << "b"
        << setw(12) << "c" << setw(16) << "dF(c)"
         << endl;
    while ((b - a) > epsilon) {
        if (FA == 0) {
            FA = MainFunction(a);
        }
        if (FB == 0) {
            FB = MainFunction(b);
        }
        if (DFA == 0) {
            DFA = MainFunctionDif(a);
        }
        if (DFB == 0) {
            DFB = MainFunctionDif(b);
        }

        if (DFA >= 0) {
            cout << "Min of Tangent method = " << setprecision(5) << a;
            break;
        }

        if (DFB <= 0) {
            cout << "Min of Tangent method = " << setprecision(5) << b;
            break;
        }

        c = (a * DFA - b * DFB + FB - FA) / (DFA - DFB);
        DFC = MainFunctionDif(c);


        cout << setw(12) << setprecision(7) << iteration << setw(16) << setprecision(7) << a
            << setw(12) << setprecision(7) << b
            << setw(12) << setprecision(7) << c << setw(16) << setprecision(7) << DFC
            << endl;
        if (DFC > 0) {
            a = a;
            b = c;
            FB = DFB = 0;
        }
        if (DFC < 0) {
            a = c;
            b = b;
            FA = DFA = 0;
        }
        if (DFC == 0) {
            cout << "Min of Tangent method = " << setprecision(5) << c;
            break;
        }

        iteration += 1;
    }
    cout << "Min of Tangent method = " << setprecision(4) << (a + b) / 2 << endl;
    cout << "Number of function calls = " << Count << endl;
    Count = 0;

    cout << "Newton - Raphson method method" << endl;
    a = -3;
    DFA = 111110;
    double DDFA;
    DDFA = 0;
    iteration = 1;
    cout << setw(12) << "Iteration " << setw(19) << "a "
        << setw(16) << "dF(a)"
        << setw(16) << "ddF(a)" 
        << endl;
    while (abs(DFA) > epsilon) {
        DFA = MainFunctionDif(a);
        DDFA = MainFunctionDif2(a);



        cout << setw(12) << setprecision(7) << iteration << setw(19) << setprecision(7) << a
            << setw(16) << setprecision(7) << DFA
            << setw(16) << setprecision(7) << DDFA
            << endl;



        a = a - DFA / DDFA;
        iteration += 1;
    }
    cout << "Min of Tangent method = " << setprecision(4) << a << endl;
    cout << "Number of function calls = " << Count << endl;
    
    cout << MainFunctionDif(-5);
    return 0;
}