#include <iostream>                              
#include <cmath>   
#include <iomanip>
#include <complex>
#include <algorithm> 

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
	double Modul = std::abs((z - C1) * (z - C2) * (z - C3));
	return Modul * Modul;
}
double FI(complex< double > z, complex< double > C1, complex< double > C2, complex< double > C3, double r) {
	
	double Modul = Func(z,C1,C2,C3) + r * (std::max(0.0,real(z)*real(z) + imag(z)*imag(z) - 16)) * (std::max(0.0, real(z) * real(z) + imag(z) * imag(z) - 16));
	return Modul * Modul;
}
double FI_2(complex< double > z, complex< double > C1, complex< double > C2, complex< double > C3, double t) {
	
	double Modul = Func(z, C1, C2, C3) - t * ( 1 / (real(z) * real(z) + imag(z) * imag(z) - 16) );
	return Modul * Modul;
}
complex<double> Gradient(complex< double > z, double r) {
	double xx, yy, h;
	/*xx = DFunc_X(real(z), imag(z));
	yy = DFunc_Y(real(z), imag(z));*/
	h = 0.000000001;
	complex<double> hx(h, 0);
	complex<double> hy(0, h);
	double Fz, Fzy, Fzx;
	Fz = FI(z, C1, C2, C3, r);
	Fzx = FI(z + hx, C1, C2, C3, r);
	Fzy = FI(z + hy, C1, C2, C3, r);
	xx = (-Fz + Fzx) / h;
	yy = (Fzy - Fz) / h;
	complex<double> Grad(xx, yy);
	//cout << "Grad" << Grad;
	Count -= 1;
	return Grad;
}
complex<double> Gradient_2(complex< double > z, double t) {
	double xx, yy, h;
	/*xx = DFunc_X(real(z), imag(z));
	yy = DFunc_Y(real(z), imag(z));*/
	h = 0.000000001;
	complex<double> hx(h, 0);
	complex<double> hy(0, h);
	double Fz, Fzy, Fzx;
	Fz = FI_2(z, C1, C2, C3, t);
	Fzx = FI_2(z + hx, C1, C2, C3, t);
	Fzy = FI_2(z + hy, C1, C2, C3, t);
	xx = (-Fz + Fzx) / h;
	yy = (Fzy - Fz) / h;
	complex<double> Grad(xx, yy);
	//cout << "Grad" << Grad;
	Count -= 1;
	return Grad;
}
double FuncMNGS(complex< double > Grad, complex< double > x, double alpha, double r){
	Count += 1;
	complex <double> Grad_N;
	Grad_N = Gradient(x - alpha * Grad, r);
	Grad_N = Grad_N / sqrt(real(Grad_N) * real(Grad_N) + imag(Grad_N) * imag(Grad_N));
	double Modul = -real(Grad) * real(Grad_N) + -imag(Grad) * imag(Grad_N);
	return Modul;
}
double FuncMNGS_2(complex< double > Grad, complex< double > x, double alpha, double t) {
	Count += 1;
	complex <double> Grad_N;
	Grad_N = Gradient_2(x - alpha * Grad, t);
	Grad_N = Grad_N / sqrt(real(Grad_N) * real(Grad_N) + imag(Grad_N) * imag(Grad_N));
	double Modul = -real(Grad) * real(Grad_N) + -imag(Grad) * imag(Grad_N);
	return Modul;
}
complex<double> Gradient_Func(complex< double > z) {
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
complex < double> Gradient_G(complex<double> z) {
	complex<double> Grad(2*real(z),2*imag(z));
	
	return Grad;
}
complex<double> Opt_X(complex< double > z, double H1, double H2, double r) {
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
		F_New = FI(z, C1, C2, C3,r);
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
		
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
			<< endl;

	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H11;
			F_New = FI(z, C1, C2, C3,r);
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
			
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << -(H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
				<< endl;
		};
	};
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}
complex<double> Opt_Y(complex< double > z, double H1, double H2,double r) {
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
		F_New = FI(z, C1, C2, C3,r);
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
		
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << (H2) << setw(18) << setprecision(8) << F_New
			<< endl;
	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H22;
			F_New = FI(z, C1, C2, C3,r);
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
			
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << -(H2) << setw(18) << setprecision(8) << F_New
				<< endl;
		};
	}
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}
complex<double> Opt_X_2(complex< double > z, double H1, double H2, double t) {
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
		F_New = FI_2(z, C1, C2, C3, t);
		if ((real(z) * real(z) + imag(z) * imag(z) - 16) < 0) {
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
		
		
		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
			<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
			<< endl;
		}
		else {
			stopper = 1;
			F_otvet = F;
			F = F_New;
			otvet = z - H11;
			break;
		}

	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H11;
			F_New = FI_2(z, C1, C2, C3, t);
			if ((real(z) * real(z) + imag(z) * imag(z) - 16) < 0) {
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
			
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << -(H1) << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << F_New
				<< endl;
			}
			else {
				otvet = z + H11;
				F_otvet = F;
			}
		};
	};
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}
complex<double> Opt_Y_2(complex< double > z, double H1, double H2, double t) {
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
		F_New = FI_2(z, C1, C2, C3, t);
		if ((real(z) * real(z) + imag(z) * imag(z) - 16) < 0) {
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

			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
				<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << (H2) << setw(18) << setprecision(8) << F_New
				<< endl;
		}
		else {
			stopper = 1;
			F_otvet = F;
			F = F_New;
			otvet = z - H22;
		}
	}
	z = start;
	F = F_start;
	if (triada == 0) {
		while (stopper > 0) {
			z = z - H22;
			F_New = FI_2(z, C1, C2, C3, t);
			if ((real(z) * real(z) + imag(z) * imag(z) - 16) < 0) {
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

				cout << setw(9) << iteration << setw(12) << setprecision(8) << real(z)
					<< setw(12) << setprecision(8) << imag(z) << setw(19) << setprecision(8) << (H1) << setw(19) << setprecision(8) << -(H2) << setw(18) << setprecision(8) << F_New
					<< endl;
			}
			else {
				otvet = z + H22;
				F_otvet = F;
			}
		};
	}
	//cout << z << endl;
	F_ret = F_otvet;
	return otvet;
}
int main(){
	complex<double> x(6,15);
	complex<double> Grad, New_x, First_root, grad, ogr;
	double r = 1;
	double globalEPS = 0.00005;
	double POL;
	double g, Fi, angle;
	Count = 0;
	std::cout << setw(9) << "Iteration" << setw(30) << setprecision(8) << "(x,y)"  << setw(50) << setprecision(8) << "Gradient F"
		<< setw(24) << "|F(x,y)|^2" << setw(24) << "g(x,y)" << setw(24) << "FI(x,y)" << setw(20) << "ANGLE"
		<< endl;
	iteration = 0;
	Grad = Gradient(x, r);
	while ((std::max(0.0, real(x) * real(x) + imag(x) * imag(x) - 4)) * (std::max(0.0, real(x) * real(x) + imag(x) * imag(x) - 16)) > globalEPS){
		
		double shag = 10;
		double alpha = 10, start_alpha = 10;
		double F, F_new, Gradient_epsilon, lambda;
		Gradient_epsilon = 0.1;
		lambda = 0.35;
		iteration = 0;
		start_alpha = alpha = 0.0009;
	
		g = real(x) * real(x) + imag(x) * imag(x) - 16;
		/*POL = Func(x, C1, C2, C3);
		g = real(x) * real(x) + imag(x) * imag(x) - 16;
		F = FI(x, C1, C2, C3, r);*/
		/*grad = Gradient_Func(x);
		ogr = Gradient_G(x);
		angle = acos(( real(grad) * real(ogr) + imag(grad) * imag(ogr) )/ sqrt(real(grad) * real(grad) + imag(grad) * imag(grad)) / sqrt(real(ogr) * real(ogr) + imag(ogr) * imag(ogr)));
		std::cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << grad
			<< setw(24) << POL << setw(24) << g << setw(24) << F << setw(20) << angle
			<< endl;*/
		Count -= 3;
		F = FI(x, C1, C2, C3, r);

		/*cout << abs(Grad) << endl << (real(Grad) * real(Grad) + imag(Grad) * imag(Grad)) << endl;*/
		
		/*cout << setw(9) << "Iteration " << setw(30) << "(x,y)" << setw(50) << "Gradient" <<
			setw(30) << "alpha" << setw(24) << "|F(x,y)|^2"
			<< endl;*/
		double epsilon = 0.0005;
		double H1, H2, F_start;
		H1 = H2 = 5;
		
		complex< double > New_min, Prev_min, prev;
		Prev_min = x;
		New_min = x;
		iteration = 0;
		F_ret = FI(x, C1, C2, C3,r);
		F_start = F_ret;
		cout << setw(9) << "Iteration " << setw(12) << "x"
			<< setw(12) << "y" << setw(12) << "H1" << setw(12) << "H2" << setw(24) << "|F(x,y)|^2"
			<< endl;
		/*cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
			<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
			<< endl;
		Count -= 1;*/
		New_min = Opt_X(Prev_min, H1, H2,r);
		H1 /= 2.1;

		New_min = Opt_Y(New_min, H1, H2,r);
		H2 /= 2.1;

		/*iteration += 1;

		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
			<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
			<< endl;
		Count -= 1;*/
		while ((/*abs(New_min - Prev_min) > epsilon) ||*/ H1 > epsilon)) {

			prev = Prev_min;
			Prev_min = New_min;
			New_min = Opt_X(Prev_min, H1, H2,r);
			H1 /= 2.1;
			New_min = Opt_Y(New_min, H1, H2,r);
			H2 /= 2.1;
			/*iteration += 1;
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
				<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
				<< endl;
			Count -= 1;*/
			//cout << setprecision(8) << Prev_min << endl;
			//cout << setprecision(8) << New_min << endl;
		};
		x = New_min;
		iteration += 1;
		/*std::cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad
			<< setw(24) << POL
			<< endl;*/
		
		
		r = 10 * r;
		POL = Func(x, C1, C2, C3);
		g = real(x) * real(x) + imag(x) * imag(x) - 16;
		grad = Gradient_Func(x);
		ogr = Gradient_G(x);
		angle = ((real(grad) * real(ogr) + imag(grad) * imag(ogr)) / sqrt(real(grad) * real(grad) + imag(grad) * imag(grad)) / sqrt(real(ogr) * real(ogr) + imag(ogr) * imag(ogr)));
		/*angle = angle / 3.14 * 180;*/
		Count -= 3;
		std::cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Gradient_Func(x)
			<< setw(24) << POL << setw(24) << g << setw(24) << F << setw(20) << angle
			<< endl;
	}
	
	
	std::cout << "Uslov. min: " << x << endl;
	std::cout << "Number of function calls: " << Count << endl;
	Count = 0;
	complex <double> xx(0, 0);
	x = xx;
	double t = 1;
	std::cout << setw(9) << "Iteration" << setw(30) << setprecision(8) << "(x,y)" << setw(50) << setprecision(8) << "Gradient F"
		<< setw(24) << "|F(x,y)|^2" << setw(24) << "g(x,y)" << setw(24) << "FI(x,y)"
		<< endl;
	
	Grad = Gradient_2(x, t);
	while ( abs(t * ( - 1 / (real(x) * real(x) + imag(x) * imag(x) - 16  ) )) > globalEPS) {
		
		double shag ;
		double alpha , start_alpha;
		alpha = start_alpha = shag = 0.35;
		double F, F_new, Gradient_epsilon, lambda;
		/*Grad = Gradient_2(x, t);
		Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));*/
		g = real(x) * real(x) + imag(x) * imag(x) - 16;
		/*F = FuncMNGS_2(Grad, x, alpha, t);*/
		/*POL = Func(x, C1, C2, C3);
		g = real(x) * real(x) + imag(x) * imag(x) - 16;
		Fi = FI_2(x, C1, C2, C3, t);
		std::cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Gradient_Func(x)
			<< setw(24) << POL << setw(24) << g << setw(24) << Fi
			<< endl;*/
		
		Gradient_epsilon = 0.1;
		lambda = 0.35;
		
		start_alpha = alpha = 0.35;
		

		/*cout << abs(Grad) << endl << (real(Grad) * real(Grad) + imag(Grad) * imag(Grad)) << endl;*/
		F = FI_2(x, C1, C2, C3, t);
		/*cout << setw(9) << "Iteration " << setw(30) << "(x,y)" << setw(50) << "Gradient" <<
			setw(30) << "alpha" << setw(24) << "|F(x,y)|^2"
			<< endl;*/
		double epsilon = 0.0005;
		double H1, H2, F_start;
		H1 = H2 = 0.5;

		complex< double > New_min, Prev_min, prev;
		Prev_min = x;
		New_min = x;
		iteration = 0;
		F_ret = FI_2(x, C1, C2, C3, r);
		F_start = F_ret;
		cout << setw(9) << "Iteration " << setw(12) << "x"
			<< setw(12) << "y" << setw(12) << "H1" << setw(12) << "H2" << setw(24) << "|F(x,y)|^2"
			<< endl;
		/*cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
			<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
			<< endl;
		Count -= 1;*/
		New_min = Opt_X_2(Prev_min, H1, H2, t);
		H1 /= 2.1;

		New_min = Opt_Y_2(New_min, H1, H2, t);
		H2 /= 2.1;

		/*iteration += 1;

		cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
			<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
			<< endl;
		Count -= 1;*/
		while (/*(abs(New_min - Prev_min) > epsilon) && */(H1 > epsilon)) {

			prev = Prev_min;
			Prev_min = New_min;
			New_min = Opt_X_2(Prev_min, H1, H2, t);
			H1 /= 2.1;
			New_min = Opt_Y_2(New_min, H1, H2, t);
			H2 /= 2.1;
			/*iteration += 1;
			cout << setw(9) << iteration << setw(12) << setprecision(8) << real(New_min)
				<< setw(12) << setprecision(8) << imag(New_min) << setw(19) << setprecision(8) << H1 << setw(19) << setprecision(8) << H2 << setw(18) << setprecision(8) << Func(New_min, C1, C2, C3)
				<< endl;
			Count -= 1;*/
			//cout << setprecision(8) << Prev_min << endl;
			//cout << setprecision(8) << New_min << endl;
		};
		x = New_min;

		
		iteration += 1;
		/*std::cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Grad
			<< setw(24) << POL
			<< endl;*/

		t = t / 5.8;
		
		POL = Func(x, C1, C2, C3);
		g = real(x) * real(x) + imag(x) * imag(x) - 16;
		Fi = FI_2(x, C1, C2, C3, t);
		Count -= 5;
		grad = Gradient_Func(x);
		ogr = Gradient_G(x);
		angle = ((real(grad) * real(ogr) + imag(grad) * imag(ogr)) / sqrt(real(grad) * real(grad) + imag(grad) * imag(grad)) / sqrt(real(ogr) * real(ogr) + imag(ogr) * imag(ogr)));
		/*angle = angle / 3.14 * 180;*/
		Count -= 1;
		std::cout << setw(9) << iteration << setw(30) << setprecision(8) << x << setw(50) << setprecision(8) << Gradient_Func(x)
			<< setw(24) << POL << setw(24) << g << setw(24) << F << setw(20) << angle
			<< endl;
	}


	std::cout << "Uslov. min: " << x << endl;
	std::cout << "Number of function calls: " << Count << endl;



	
	return 0;
}