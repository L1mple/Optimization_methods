#include <iostream>                              
#include <cmath>   
#include <iomanip>
#include <complex>

using namespace std;

static double Name(4.0);
static double Surname(5.0);
static double Patronymic(10.0);
static double Count(0);
static double iteration(0);
static double F_ret(0);
complex< double > C1(Name, Patronymic);
complex< double > C2(Patronymic, Name);
complex< double > C3(Name, Surname);

//double Func(complex< double > z) {
//	Count += 1;
//	double Modul = std::abs(z*z*z - (C1 + C2 + C3)*z*z + z*(C1 * C2 + C2 * C3 + C3 * C1) - (C1 * C2 * C3));
//	return Modul * Modul;
//}
double Func(complex< double > z, complex< double > C1, complex< double > C2, complex< double > C3) {
	Count += 1;
	double Modul = std::abs((z - C1)*(z - C2)*(z - C3));
	return Modul * Modul;
}
double Func_2(complex< double > z, complex< double > C1, complex< double > C2, complex< double > R) {
	Count += 1;
	double Modul = std::abs(z*z + C1*z + C2 + R);
	return Modul * Modul;
}
double Func_3(complex< double > z, complex< double > C1, complex< double > R) {
	Count += 1;
	double Modul = std::abs((z + C1 + R));
	return Modul * Modul;
}
double DFunc_X(double x, double y) {
	double Modul = 2*(x*x*x - 3*y*y*x + 18*y*y - 18*x*x + 38*x*y - 14*x - 242*y - 580)*(3*x*x - 3*y*y - 36*x + 38*y - 14) + 2*(3*x*x*y - y*y*y - 36*x*y - 19*x*x + 19*y*y + 242*x - 14*y + 464)*(3*2*x*y - 36*y - 38*x + 242);
	return Modul;
}
double DFunc_Y(double x, double y) {
	double Modul = 2*(x * x * x - 3 * y * y * x + 18 * y * y - 18 * x * x + 38 * x * y - 14 * x - 242 * y - 580)*(-6*y*x + 36*y + 38*x - 242) + 2*(3 * x * x * y - y * y * y - 36 * x * y - 19 * x * x + 19 * y * y + 242 * x - 14 * y + 464)*(3*x*x - 3*y*y - 36*x + 38*y - 14);
		return Modul;
}
//complex<double> Gradient(complex< double > z) {
//	double xx, yy;
//	xx = DFunc_X(real(z), imag(z));
//	yy = DFunc_Y(real(z), imag(z));
//	
//	complex<double> Grad(xx, yy);
//	//cout << "Grad" << Grad;
//	Count += 1;
//	return Grad;
//}
complex<double> Gradient(complex< double > z) {
	double xx, yy, h;
	/*xx = DFunc_X(real(z), imag(z));
	yy = DFunc_Y(real(z), imag(z));*/
	h = 0.000000001;
	complex<double> hx(h, 0);
	complex<double> hy(0, h);
	double Fz, Fzy, Fzx;
	Fz = Func(z, C1, C2, C3);
	Fzx = Func(z + hx, C1, C2, C3);
	Fzy = Func(z + hy, C1, C2, C3);
	xx = (-Fz + Fzx) / h;
	yy = (Fzy - Fz) / h;
	complex<double> Grad(xx, yy);
	//cout << "Grad" << Grad;
	Count -= 1;
	return Grad;
}
complex<double> Gradient_2(complex< double > z, complex< double > B1, complex< double > B2, complex< double > R_1) {
	double xx, yy, h;
	/*xx = DFunc_X(real(z), imag(z));
	yy = DFunc_Y(real(z), imag(z));*/
	h = 0.000000001;
	complex<double> hx(h, 0);
	complex<double> hy(0, h);
	double Fz, Fzy, Fzx;
	Fz = Func_2(z, B1, B2, R_1);
	Fzx = Func_2(z + hx, B1, B2, R_1);
	Fzy = Func_2(z + hy, B1, B2, R_1);
	xx = (-Fz + Fzx) / h;
	yy = (Fzy - Fz) / h;
	
	complex<double> Grad(xx, yy);
	//cout << "Grad" << Grad;
	Count -= 1;
	return Grad;
}
complex<double> Gradient_3(complex< double > z, complex< double > A1, complex< double > R_2) {
	double xx, yy, h;
	/*xx = DFunc_X(real(z), imag(z));
	yy = DFunc_Y(real(z), imag(z));*/
	h = 0.000000001;
	complex<double> hx(h, 0);
	complex<double> hy(0, h);
	double Fz, Fzy, Fzx;
	Fz = Func_3(z, A1, R_2);
	Fzx = Func_3(z + hx, A1, R_2);
	Fzy = Func_3(z + hy, A1, R_2);
	xx = (-Fz + Fzx) / h;
	yy = (Fzy - Fz) / h;
	
	complex<double> Grad(xx, yy);
	//cout << "Grad" << Grad;
	Count -= 1;
	return Grad;
}





complex<double> Opt_X(complex< double > z, double H1, double H2) {
	double stopper = -1;
	complex< double > H11(H1, 0);
	double F_New, F_start, F_otvet;
	complex< double > start, otvet;
	start = z;
	double F = F_ret;
	double triada;
	F_start = F;
	triada = 0;
	cout << setw(9) << "Start" << setw(12) << setprecision(8) << real(z)
		<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F
		<< endl;
	while (stopper < 0) {
		z = z + H11;
		F_New = Func(z, C1, C2, C3);
		if (F < F_New) {
			stopper = 1;
			F_otvet = F;
			F = F_New;
			otvet = z - H11;
		};
		if (F > F_New) {
			F = F_New;
			triada += 1;
		};
		iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1)<< setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
			<< endl;
		
	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H11;
			F_New = Func(z, C1, C2, C3);
			if (F < F_New) {
				stopper = 0;
				if (F < F_otvet) {
					otvet = z + H11;
					F_otvet = F;
				}
				//z = z + H11;
			};
			if (F > F_New) {
				F = F_New;
			};
			iteration += 1;
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << -(H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
				<< endl;
		};
	};
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}



complex<double> Opt_X_2(complex< double > z, double H1, double H2, complex< double > B1, complex< double > B2, complex< double > R) {
	double stopper = -1;
	complex< double > H11(H1, 0);
	double F_New, F_start, F_otvet;
	complex< double > start, otvet;
	start = z;
	double F = F_ret;
	double triada;
	F_start = F;
	triada = 0;
	cout << setw(9) << "Start" << setw(12) << setprecision(8) << real(z)
		<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F
		<< endl;
	while (stopper < 0) {
		z = z + H11;
		F_New = Func_2(z, B1, B2, R);
		if (F < F_New) {
			stopper = 1;
			F_otvet = F;
			F = F_New;
			otvet = z - H11;
		};
		if (F > F_New) {
			F = F_New;
			triada += 1;
		};
		iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
			<< endl;
	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H11;
			F_New = Func_2(z, B1, B2, R);
			if (F < F_New) {
				stopper = 0;
				if (F < F_otvet) {
					otvet = z + H11;
					F_otvet = F;
				}
				//z = z + H11;
			};
			if (F > F_New) {
				F = F_New;
			};
			iteration += 1;
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << -(H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
				<< endl;
		};
	}
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}

complex<double> Opt_X_3(complex< double > z, double H1, double H2, complex< double > A1, complex< double > R) {
	double stopper = -1;
	complex< double > H11(H1, 0);
	double F_New, F_start, F_otvet;
	complex< double > start, otvet;
	start = z;
	double F = F_ret;
	double triada;
	F_start = F;
	triada = 0;
	cout << setw(9) << "Start" << setw(12) << setprecision(8) << real(z)
		<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F
		<< endl;
	while (stopper < 0) {
		z = z + H11;
		F_New = Func_3(z, A1, R);
		if (F < F_New) {
			stopper = 1;
			F_otvet = F;
			F = F_New;
			otvet = z - H11;
		};
		if (F > F_New) {
			F = F_New;
			triada += 1;
		};
		iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
			<< endl;
	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H11;
			F_New = Func_3(z, A1, R);
			if (F < F_New) {
				stopper = 0;
				if (F < F_otvet) {
					otvet = z + H11;
					F_otvet = F;
				}
				//z = z + H11;
			};
			if (F > F_New) {
				F = F_New;
			};
			iteration += 1;
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << -(H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
				<< endl;
		};
	}
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}







complex<double> Opt_Y(complex< double > z, double H1, double H2) {
	double stopper = -1;
	complex< double > H22(0, H2);
	double F_New, F_start, F_otvet;
	complex< double > start, otvet;
	start = z;
	double F = F_ret;
	double triada;
	F_start = F;
	triada = 0;
	cout << setw(9) << "Start" << setw(12) << setprecision(8) << real(z)
		<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F
		<< endl;
	while (stopper < 0) {
		z = z + H22;
		F_New = Func(z, C1, C2, C3);
		if (F < F_New) {
			stopper = 1;
			F_otvet = F;
			F = F_New;
			otvet = z - H22;
		};
		if (F > F_New) {
			F = F_New;
			triada += 1;
		};
		iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << (H2) << setw(18) << setprecision(8) << F_New
			<< endl;
	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H22;
			F_New = Func(z, C1, C2, C3);
			if (F < F_New) {
				stopper = 0;
				if (F < F_otvet) {
					otvet = z + H22;
					F_otvet = F;
				};
				//z = z + H22;
			};
			if (F > F_New) {
				F = F_New;
			};
			iteration += 1;
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << -(H2) << setw(18) << setprecision(8) << F_New
				<< endl;
		};
	}
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}


complex<double> Opt_Y_2(complex< double > z, double H1, double H2, complex< double > B1, complex< double > B2, complex< double > R) {
	double stopper = -1;
	complex< double > H22(0, H2);
	double F_New, F_start, F_otvet;
	complex< double > start, otvet;
	start = z;
	double F = F_ret;
	double triada;
	F_start = F;
	triada = 0;
	cout << setw(9) << "Start" << setw(12) << setprecision(8) << real(z)
		<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F
		<< endl;
	while (stopper < 0) {
		z = z + H22;
		F_New = Func_2(z, B1, B2, R);
		if (F < F_New) {
			stopper = 1;
			F_otvet = F;
			F = F_New;
			otvet = z - H22;
		};
		if (F > F_New) {
			F = F_New;
			triada += 1;
		};
		iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << (H2) << setw(18) << setprecision(8) << F_New
			<< endl;
	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H22;
			F_New = Func_2(z, B1, B2, R);
			if (F < F_New) {
				stopper = 0;
				if (F < F_otvet) {
					otvet = z + H22;
					F_otvet = F;
				}
				//z = z + H22;
			};
			if (F > F_New) {
				F = F_New;
			};
			iteration += 1;
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << -(H2) << setw(18) << setprecision(8) << F_New
				<< endl;
		};
	}
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}

complex<double> Opt_Y_3(complex< double > z, double H1, double H2, complex< double > A1, complex< double > R) {
	double stopper = -1;
	complex< double > H22(0, H2);
	double F_New, F_start, F_otvet;
	complex< double > start, otvet;
	start = z;
	double F = F_ret;
	double triada;
	F_start = F;
	triada = 0;
	cout << setw(9) << "Start" << setw(12) << setprecision(8) << real(z)
		<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F
		<< endl;
	while (stopper < 0) {
		z = z + H22;
		F_New = Func_3(z, A1, R);
		if (F < F_New) {
			stopper = 1;
			F_otvet = F;
			F = F_New;
			otvet = z - H22;
		};
		if (F > F_New) {
			F = F_New;
			triada += 1;
		};
		iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << (H2) << setw(18) << setprecision(8) << F_New
			<< endl;
	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H22;
			F_New = Func_3(z, A1, R);
			if (F < F_New) {
				stopper = 0;
				if (F < F_otvet) {
					otvet = z + H22;
					F_otvet = F;
				}
				//z = z + H22;
			};
			if (F > F_New) {
				F = F_New;
			};
			iteration += 1;
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << -(H2) << setw(18) << setprecision(8) << F_New
				<< endl;
		};
	}
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}


double FuncMNGS(complex< double > Grad, complex< double > x, double alpha) {
	
	complex <double> Grad_N;
	Grad_N = Gradient(x - alpha * Grad);
	Grad_N = Grad_N / sqrt(real(Grad_N) * real(Grad_N) + imag(Grad_N) * imag(Grad_N));
	double Modul = abs(-real(Grad) * real(Grad_N) + -imag(Grad) * imag(Grad_N));
	return Modul;
}

double FuncMNGS_2(complex< double > Grad, complex< double > x, double alpha, complex< double > B1, complex< double > B2, complex< double > R_1) {
	Count += 1;
	complex <double> Grad_N;
	Grad_N = Gradient_2(x - alpha * Grad, B1, B2, R_1);
	Grad_N = Grad_N / sqrt(real(Grad_N) * real(Grad_N) + imag(Grad_N) * imag(Grad_N));
	double Modul = -real(Grad) * real(Grad_N) + -imag(Grad) * imag(Grad_N);
	return Modul;
}


double FuncMNGS_3(complex< double > Grad, complex< double > x, double alpha, complex< double > A1, complex< double > R_2) {
	Count += 1;
	complex <double> Grad_N;
	Grad_N = Gradient_3(x - alpha * Grad, A1, R_2);
	Grad_N = Grad_N / sqrt(real(Grad_N) * real(Grad_N) + imag(Grad_N) * imag(Grad_N));
	double Modul = -real(Grad) * real(Grad_N) + -imag(Grad) * imag(Grad_N);
	return Modul;
}






double GR(complex<double> x,complex<double> Grad, double alpha)//
{
	double a = 0.0000001;
	double b = alpha;
	double s5 = sqrt(5);
	double c = (3 - s5) / 2 * (b - a) + a;
	double d = (s5 - 1) / 2 * (b - a) + a;
	while ((b - a) / 2 > 0.00005)
	{
		if (FuncMNGS(Grad, x, c) <= FuncMNGS(Grad, x, d))
		{
			b = d;
			d = c;
			c = (3 - s5) / 2 * (b - a) + a;
		}
		else
		{
			a = c;
			c = d;
			d = (s5 - 1) / 2 * (b - a) + a;
		}
	}
	return (a + b) / 2;
}








int main() {
	complex<double> First_root, Second_root, Third_root;
	double epsilon = 0.00005;	
	double H1, H2, F_start;
	H1 = H2 = 5;
	complex< double > minimum(0, 0);
	complex< double > New_min, Prev_min, prev;
	Prev_min = minimum;
	New_min = minimum;
	iteration = 0; 
	F_ret = Func(minimum, C1, C2, C3);
	F_start = F_ret;
	cout << setw(9) << "Iteration " << setw(12) << "x"
		<< setw(12) << "y" << setw(12) << "H1" << setw(12) << "H2" << setw(24) << "|F(x,y)|^2"
		<< endl;
	/*cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
		<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
		<< endl;
	Count -= 1;*/
	New_min = Opt_X(Prev_min, H1, H2);
	H1 /= 2.1;
	
	New_min = Opt_Y(New_min, H1, H2);
	H2 /= 2.1;
	
	/*iteration += 1;
	
	cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
		<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
		<< endl;
	Count -= 1;*/
	while ((abs(New_min - Prev_min) > epsilon) | (H1 > epsilon)) {		
		
		prev = Prev_min;
		Prev_min = New_min;
		New_min = Opt_X(Prev_min, H1, H2);
		H1 /= 2.1;
		New_min = Opt_Y(New_min, H1, H2);
		H2 /= 2.1;
		/*iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
			<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
			<< endl;
		Count -= 1;*/
		//cout << setprecision(8) << Prev_min << endl;
		//cout << setprecision(8) << New_min << endl;
	};
	First_root = New_min;
	cout << "First root: " << setprecision(8) << New_min << endl;
	// 	Horner sheme
	complex< double > B1, B2, R_1;
	B1 = New_min + -C1 - C2 - C3;
	B2 = New_min * B1 + (C1 * C2 + C2 * C3 + C3 * C1);
	R_1 = B2 * New_min + (-C1 * C2 * C3);
	cout << "Function after Horner: |z^2 + z * " << B1 << " + " << B2 << "+ " << R_1 << "|^2" << endl << endl;
	//
	H1 = H2 = 5;
	minimum = (0, 0);
	Prev_min = minimum;
	New_min = minimum;
	iteration = 0;
	F_ret = F_start;
	cout << setw(9) << "Iteration " << setw(12) << "x"
		<< setw(12) << "y" << setw(12) << "H1" << setw(12) << "H2" << setw(24) << "|F(x,y)|^2"
		<< endl;
	/*cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
		<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
		<< endl;
	Count -= 1;*/
	New_min = Opt_X_2(Prev_min, H1, H2, B1, B2, R_1);
	H1 /= 2.1;
	New_min = Opt_Y_2(New_min, H1, H2, B1, B2, R_1);
	H2 /= 2.1;
	/*iteration += 1;
	cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
		<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
		<< endl;
	Count -= 1;*/
	while ((abs(New_min - Prev_min) > epsilon) | (H1 > epsilon)) {
		prev = Prev_min;
		Prev_min = New_min;
		New_min = Opt_X_2(New_min, H1, H2, B1, B2, R_1);
		H1 /= 2.1;
		New_min = Opt_Y_2(New_min, H1, H2, B1, B2, R_1);
		H2 /= 2.1;
		/*iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
			<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
			<< endl;
		Count -= 1;*/
		//cout << setprecision(8) << Prev_min << endl;
		//cout << setprecision(8) << New_min << endl;
	};
	Second_root = New_min;
	cout << "Second root: " << New_min << endl;

	// Horner sheme
	complex< double > A1, R_2;
	A1 = New_min + B1;
	R_2 = New_min * A1 + B2 + R_1;
	cout << "Function after Horner: |z - " << A1 << "  + " << R_2 << "|^2" << endl;
	//
	H1 = H2 = 5;
	minimum = (0, 0);
	Prev_min = minimum;
	New_min = minimum;
	iteration = 0;
	F_ret = F_start;
	cout << setw(9) << "Iteration " << setw(12) << "x"
		<< setw(12) << "y" << setw(12) << "H1" << setw(12) << "H2" << setw(24) << "|F(x,y)|^2"
		<< endl;
	//cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
	//	<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
	//	<< endl;
	//Count -= 1;
	New_min = Opt_X_3(Prev_min, H1, H2, A1, R_2);
	H1 /= 2.1;
	New_min = Opt_Y_3(New_min, H1, H2, A1, R_2);
	H2 /= 2.1;
	/*iteration += 1;
	cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
		<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
		<< endl;
	Count -= 1;*/
	while ((abs(New_min - Prev_min) > epsilon) | (H1 > epsilon)) {
		prev = Prev_min;
		Prev_min = New_min;
		New_min = Opt_X_3(New_min, H1, H2, A1, R_2);
		H1 /= 2.1;
		New_min = Opt_Y_3(New_min, H1, H2, A1, R_2);
		H2 /= 2.1;
		/*iteration += 1;
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
			<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
			<< endl;
		Count -= 1;*/
		//cout << setprecision(8) << Prev_min << endl;
		//cout << setprecision(8) << New_min << endl;
	};
	Third_root = New_min;
	cout << "Third root: " << New_min << endl << endl;
	cout << "Actual roots: " << endl << C1 << endl << C2 << endl << C3 << endl << endl;
	cout << "Roots from function :" << endl;
	cout << "First root = " << First_root << endl;
	cout << "Second root = " << Second_root << endl;
	cout << "Third root = " << Third_root << endl << endl;
	cout << "Number of function calls: " << Count;

	cout << endl << endl << endl << endl;
	cout << "Gradient" << endl;
	cout << "Pitch splitting" << endl;
	Count = 0;
	double Gradient_epsilon, lambda, alpha, start_alpha, F_new, F;
	complex <double> x(4.5, 4.5);
	complex <double> New_x(4.5, 4.5);
	complex <double> Grad;
	
	Gradient_epsilon = 0.1;
	lambda = 0.35;
	iteration = 0;
	start_alpha = alpha  = 0.0009;
	Grad = Gradient(x);
	
	cout << abs(Grad) << endl << (real(Grad) * real(Grad) + imag(Grad) * imag(Grad)) << endl;
	F = Func(x, C1, C2, C3);
	cout << setw(9) << "Iteration " << setw(30) << "(x,y)" << setw(50) << "Gradient" <<
		 setw(30) << "alpha" << setw(24) << "|F(x,y)|^2" 
		<< endl;
	while (sqrt((real(Grad) * real(Grad) + imag(Grad) * imag(Grad))) > 0.0005) {
		New_x = x - alpha * Grad;
		F_new = Func(New_x, C1, C2, C3);
		if ( (F_new - F) > (- alpha * Gradient_epsilon * ( real(Grad) * real(Grad) + imag(Grad) * imag(Grad) ) ) ) {
			alpha *= lambda;
			continue;
			New_x = x - alpha * Grad;
			F_new = Func(New_x, C1, C2, C3);
			
			/*cout << alpha << endl;
			cout << New_x << endl;*/
		}
		iteration += 1;
		cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
			setw(30) << alpha << setw(24) << F 
			<< endl;
		x = New_x;
		New_x = x;
		F = F_new;
		alpha = start_alpha;
		Grad = Gradient(New_x);
		/*cout << "asdasdasd" << abs(Gradient(New_x)) << endl;*/

	}
	iteration += 1;
	cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
		setw(30) << alpha << setw(24) << F
		<< endl;
	New_min = x;
	First_root = x;
	cout << "First root: " << First_root << endl;
	cout << "Number of function calls: " << Count << endl;
	Count = 0;
	// Horner
	
	B1 = New_min + -C1 - C2 - C3;
	B2 = New_min * B1 + (C1 * C2 + C2 * C3 + C3 * C1);
	R_1 = B2 * New_min + (-C1 * C2 * C3);
	cout << "Function after Horner: |z^2 + z * " << B1 << " + " << B2 << "+ " << R_1 << "|^2" << endl << endl;
	complex <double> xx(6, 20);
	//x = xx;
	//iteration = 0;
	//start_alpha = alpha = 0.009;
	//Grad = Gradient_2(x, B1, B2, R_1);
	//cout << abs(Grad) << endl << (real(Grad) * real(Grad) + imag(Grad) * imag(Grad)) << endl;
	//F = Func_2(x, B1, B2, R_1);
	//cout << setw(9) << "Iteration " << setw(30) << "(x,y)" << setw(50) << "Gradient" <<
	//	setw(30) << "alpha" << setw(24) << "|F(x,y)|^2"
	//	<< endl;
	//while (sqrt((real(Grad) * real(Grad) + imag(Grad) * imag(Grad))) > 0.0005) {
	//	New_x = x - alpha * Grad;
	//	F_new = Func_2(New_x, B1, B2, R_1);
	//	if ((F_new - F) > (+alpha * Gradient_epsilon * (real(Grad) * real(Grad) + imag(Grad) * imag(Grad)))) {
	//		alpha *= lambda;
	//		continue;
	//		New_x = x - alpha * Grad;
	//		F_new = Func_2(New_x, B1, B2, R_1);

	//		/*cout << alpha << endl;
	//		cout << New_x << endl;*/
	//	}
	//	iteration += 1;
	//	cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
	//		setw(30) << alpha << setw(24) << F
	//		<< endl;
	//	x = New_x;
	//	New_x = x;
	//	F = F_new;
	//	alpha = start_alpha;
	//	Grad = Gradient_2(New_x, B1, B2, R_1);
	//	/*cout << "asdasdasd" << abs(Gradient(New_x)) << endl;*/

	//}
	//iteration += 1;
	//cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
	//	setw(30) << alpha << setw(24) << F
	//	<< endl;
	//New_min = x;
	//Second_root = x;
	//cout << "Second root: " << Second_root << endl;
	//cout << "Number of function calls: " << Count << endl;
	//Count = 0;
	////Horner scheme
	//A1 = New_min + B1;
	//R_2 = New_min * A1 + B2 + R_1;
	//cout << "Function after Horner: |z - " << A1 << "  + " << R_2 << "|^2" << endl;
	complex <double> xxx(15, 8);
	//x = xxx;
	//iteration = 0;
	//start_alpha = alpha = 0.09;
	//Grad = Gradient_3(x, A1, R_2);
	//cout << abs(Grad) << endl << (real(Grad) * real(Grad) + imag(Grad) * imag(Grad)) << endl;
	//F = Func_3(x, A1, R_2);
	//cout << setw(9) << "Iteration " << setw(30) << "(x,y)" << setw(50) << "Gradient" <<
	//	setw(30) << "alpha" << setw(24) << "|F(x,y)|^2"
	//	<< endl;
	//while (sqrt((real(Grad) * real(Grad) + imag(Grad) * imag(Grad))) > 0.0005) {
	//	New_x = x - alpha * Grad;
	//	F_new = Func_3(New_x, A1, R_2);
	//	if ((F_new - F) > (+alpha * Gradient_epsilon * (real(Grad) * real(Grad) + imag(Grad) * imag(Grad)))) {
	//		alpha *= lambda;
	//		continue;
	//		New_x = x - alpha * Grad;
	//		F_new = Func_3(New_x, A1, R_2);

	//		/*cout << alpha << endl;
	//		cout << New_x << endl;*/
	//	}
	//	iteration += 1;
	//	cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
	//		setw(30) << alpha << setw(24) << F
	//		<< endl;
	//	x = New_x;
	//	New_x = x;
	//	F = F_new;
	//	alpha = start_alpha;
	//	Grad = Gradient_3(New_x, A1, R_2);
	//	/*cout << "asdasdasd" << abs(Gradient(New_x)) << endl;*/

	//}
	//iteration += 1;
	//cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
	//	setw(30) << alpha << setw(24) << F
	//	<< endl;
	//New_min = x;
	//Third_root = x;
	//cout << "Third root: " << Third_root << endl<< endl;
	//cout << "Number of function calls: " << Count << endl;
	//cout << "Constant pitch" << endl;
	//complex<double> z(10, 10);
	//x = z;
	F = Func(x, C1, C2, C3);
	Grad = Gradient(x);
	iteration = 0;
	cout << setw(9) << "Iteration" << setw(30) << setprecision(8) << "x" << setw(50) << setprecision(8) << "Gradient" <<
		setw(30) << "Shag" << setw(24) << "|F(x,y)|^2"
		<< endl;
	double shag = 0.0001;

	while (sqrt((real(Grad) * real(Grad) + imag(Grad) * imag(Grad))) > 0.0005) {
		New_x = x - shag * Grad;
		F_new = Func(New_x, C1, C2, C3);
		
		
		cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
			setw(30) << shag << setw(24) << F
			<< endl;
		iteration += 1;
		x = New_x;
		New_x = x;
		F = F_new;
		alpha = start_alpha;
		Grad = Gradient(New_x);
		/*cout << "asdasdasd" << abs(Gradient(New_x)) << endl;*/

	}
	iteration += 1;
	cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
		setw(30) << shag << setw(24) << F
		<< endl;
	New_min = x;
	First_root = x;
	cout << "First root: " << x << endl;
	cout << "Number of function calls: " << Count << endl;
	Count = 0;
	//// Horner

	//B1 = New_min + -C1 - C2 - C3;
	//B2 = New_min * B1 + (C1 * C2 + C2 * C3 + C3 * C1);
	//R_1 = B2 * New_min + (-C1 * C2 * C3);
	//cout << "Function after Horner: |z^2 + z * " << B1 << " + " << B2 << "+ " << R_1 << "|^2" << endl << endl;
	//x = xx;
	//F = Func_2(x, B1, B2, R_1);
	//Grad = Gradient_2(x, B1, B2, R_1);
	//iteration = 0;
	//cout << setw(9) << "Iteration" << setw(30) << setprecision(8) << "x" << setw(50) << setprecision(8) << "Gradient" <<
	//	setw(30) << "Shag" << setw(24) << "|F(x,y)|^2"
	//	<< endl;
	//shag = 0.0001;
	//while (sqrt((real(Grad) * real(Grad) + imag(Grad) * imag(Grad))) > 0.0005) {
	//	New_x = x - shag * Grad;
	//	F_new = Func_2(x, B1, B2, R_1);

	//	iteration += 1;
	//	cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
	//		setw(30) << shag << setw(24) << F
	//		<< endl;
	//	x = New_x;
	//	New_x = x;
	//	F = F_new;
	//	alpha = start_alpha;
	//	Grad = Gradient_2(New_x, B1, B2, R_1);
	//	/*cout << "asdasdasd" << abs(Gradient(New_x)) << endl;*/

	//}
	//iteration += 1;
	//cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
	//	setw(30) << shag << setw(24) << F
	//	<< endl;
	//New_min = x;
	//Second_root = x;
	//cout << "Second root: " << Second_root << endl;
	//cout << "Number of function calls: " << Count << endl;
	//Count = 0;
	////Horner scheme
	//A1 = New_min + B1;
	//R_2 = New_min * A1 + B2 + R_1;
	//cout << "Function after Horner: |z - " << A1 << "  + " << R_2 << "|^2" << endl;
	complex<double> xxxx(20, 4);
	//x = xxxx;
	//F = Func_3(x, A1, R_2);
	//Grad = Gradient_3(x, A1, R_2);
	//iteration = 0;
	//cout << setw(9) << "Iteration" << setw(30) << setprecision(8) << "x" << setw(50) << setprecision(8) << "Gradient" <<
	//	setw(30) << "Shag" << setw(24) << "|F(x,y)|^2"
	//	<< endl;
	//shag = 0.001;
	//while (sqrt((real(Grad) * real(Grad) + imag(Grad) * imag(Grad))) > 0.0005) {
	//	New_x = x - shag * Grad;
	//	F_new = Func_3(x, A1, R_2);

	//	iteration += 1;
	//	cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
	//		setw(30) << shag << setw(24) << F
	//		<< endl;
	//	x = New_x;
	//	New_x = x;
	//	F = F_new;
	//	alpha = start_alpha;
	//	Grad = Gradient_3(New_x, A1, R_2);
	//	/*cout << "asdasdasd" << abs(Gradient(New_x)) << endl;*/

	//}
	//iteration += 1;
	//cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
	//	setw(30) << shag << setw(24) << F
	//	<< endl;
	//New_min = x;
	//Third_root = x;
	//cout << "Third root: " << Third_root << endl;
	//cout << "Number of function calls: " << Count << endl;
	//Count = 0;
	complex<double> GradNew(0,0);
	
	cout << "MNGS" << endl;
	x = xx;
	shag = 10;
	alpha = start_alpha = 10;
	Grad = Gradient(x);
	double norm = sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	iteration = 0;
	F = Func(x,C1,C2,C3);
	Count -= 1;
	double POL = 31203133;
	cout << setw(9) << "Iteration" << setw(30) << setprecision(8) << "x" << setw(50) << setprecision(8) << "Gradient" <<
		setw(30) << "Shag" << setw(24) << "|F(x,y)|^2" << setw(24) << "Calls" << setw(24) << "cos(grad)" << setw(24) << "Func_MNGS"
		<< endl;
	while (shag > 0.005) {
		
		/*alpha = GR(x, Grad);*/
		New_x = x - alpha * Grad;
		F_new = Func(New_x,C1,C2,C3);
		while (F_new > F || real(Grad)*real(GradNew) + imag(Grad)*imag(GradNew) > 0) {
			alpha *= 0.9;
			New_x = x - alpha * Grad;
			F_new = Func(New_x,C1,C2,C3);
			GradNew = Gradient(New_x);
			GradNew = GradNew / sqrt(real(GradNew) * real(GradNew) + imag(GradNew) * imag(GradNew));
		}
		/*alpha = GR(x, Grad, alpha);*/
		shag = alpha;
		POL = Func(x, C1, C2, C3);
		New_min = Gradient(New_x); 
		New_min = New_min / sqrt(real(New_min) * real(New_min) + imag(New_min) * imag(New_min));
		double cosinus = real(New_min) * real(Grad) + imag(New_min) * imag(Grad);
		cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
			setw(30) << setprecision(8) << shag << setw(24) << POL << setw(24) << Count << setw(24) << cosinus << setw(24) << F
			<< endl;
		iteration += 1;
		x = x - shag * Grad;
		Grad = Gradient(x);
		norm = sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
		Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
		F = F_new;
		/*cout << "smt" << endl;*/
		
	}
	New_min = x;
	First_root = x;
	Count -= 2;
	/*cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
		setw(30) << shag << setw(24) << POL << setw(24) << Count << setw(24) << "cos(grad)"
		<< endl;*/
	cout << "First root: " << First_root << endl;
	cout << "Number of function calls: " << Count << endl;
	Count = 0;
	// Horner

	B1 = New_min + -C1 - C2 - C3;
	B2 = New_min * B1 + (C1 * C2 + C2 * C3 + C3 * C1);
	R_1 = B2 * New_min + (-C1 * C2 * C3);
	cout << "Function after Horner: |z^2 + z * " << B1 << " + " << B2 << "+ " << R_1 << "|^2" << endl << endl;
	x = xxx;
	shag = 10;
	alpha = start_alpha = 10;
	Grad = Gradient_2(x, B1, B2, R_1);
	Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	iteration = 0;
	F = FuncMNGS_2(Grad, x, alpha, B1, B2, R_1);
	
	cout << setw(9) << "Iteration" << setw(30) << setprecision(8) << "x" << setw(50) << setprecision(8) << "Gradient" <<
		setw(30) << "Shag" << setw(24) << "|F(x,y)|^2" << setw(24) << "Calls" << setw(24) << "cos(grad)" << setw(24) << "Func_MNGS"
		<< endl;
	while (shag > 0.00005) {
		New_x = x - alpha * Grad;
		F_new = FuncMNGS_2(Grad, x, alpha, B1, B2, R_1);
		while (F_new > F) {
			alpha *= 0.95;
			New_x = x - alpha * Grad;
			F_new = FuncMNGS_2(Grad, x, alpha, B1, B2, R_1);
		}
		shag = alpha;
		POL = Func_2(x, B1, B2, R_1);
		New_min = Gradient_2(New_x, B1, B2, R_1);
		New_min = New_min / sqrt(real(New_min) * real(New_min) + imag(New_min) * imag(New_min));
		double cosinus = real(New_min) * real(Grad) + imag(New_min) * imag(Grad);
		cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
			setw(30) << setprecision(8) << shag << setw(24) << POL << setw(24) << Count << setw(24) << cosinus << setw(24) << F
			<< endl;
		iteration += 1;
		x = x - shag * Grad;
		Grad = Gradient_2(x, B1, B2, R_1);
		Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
		F = F_new;
		
	}
	New_min = x;
	Second_root = x;

	/*cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
		setw(30) << shag << setw(24) << POL << setw(24) << Count << setw(24) << "cos(grad)"
		<< endl;*/
	cout << "Second root: " << Second_root << endl;
	cout << "Number of function calls: " << Count << endl;

	Count = 0;
	//Horner scheme
	A1 = New_min + B1;
	R_2 = New_min * A1 + B2 + R_1;
	cout << "Function after Horner: |z - " << A1 << "  + " << R_2 << "|^2" << endl;
	x = xxxx;
	shag = 10;
	alpha = start_alpha = 10;
	Grad = Gradient_3(x, A1, R_2);
	Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	iteration = 0;
	F = FuncMNGS_3(Grad, x, alpha, A1, R_2);

	cout << setw(9) << "Iteration" << setw(30) << setprecision(8) << "x" << setw(50) << setprecision(8) << "Gradient" <<
		setw(30) << "Shag" << setw(24) << "|F(x,y)|^2" << setw(24) << "Calls" << setw(24) << "cos(grad)" << setw(24) << "Func_MNGS"
		<< endl;
	while (shag > 0.00005) {
		New_x = x - alpha * Grad;
		F_new = FuncMNGS_3(Grad, x, alpha, A1, R_2);
		while (F_new > F) {
			alpha *= 0.9;
			New_x = x - alpha * Grad;
			F_new = FuncMNGS_3(Grad, x, alpha, A1, R_2);
		}
		shag = alpha;
		POL = Func_3(x, A1, R_2);
		New_min = Gradient_3(New_x, A1, R_2);
		New_min = New_min / sqrt(real(New_min) * real(New_min) + imag(New_min) * imag(New_min));
		double cosinus = real(New_min) * real(Grad) + imag(New_min) * imag(Grad);
		cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
			setw(30) << setprecision(8) << shag << setw(24) << POL << setw(24) << Count << setw(24) << cosinus << setw(24) << F
			<< endl;
		iteration += 1;
		x = x - shag * Grad;
		Grad = Gradient_3(x, A1, R_2);
		Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
		F = F_new;
		alpha = start_alpha;
	}
	New_min = x;
	Third_root = x;

	/*cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad <<
		setw(30) << shag << setw(24) << POL << setw(24) << Count << setw(24) << "cos(grad)"
		<< endl;*/
	cout << "Third root: " << Third_root << endl;
	cout << "Number of function calls: " << Count << endl;

	
	
	return 0;
}


