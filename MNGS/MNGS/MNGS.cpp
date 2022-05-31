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
static double EPS(0.001);
complex< double > C1(Name, Patronymic);
complex< double > C2(Patronymic, Name);
complex< double > C3(Name, Surname);



double Func(complex< double > z, complex< double > C1, complex< double > C2, complex< double > C3) {
	Count += 1;
	double Modul = std::abs((z - C1) * (z - C2) * (z - C3));
	return Modul * Modul;
}
double Func_2(complex< double > z, complex< double > C1, complex< double > C2, complex< double > R) {
	Count += 1;
	double Modul = std::abs(z * z + C1 * z + C2 + R);
	return Modul * Modul;
}
double Func_3(complex< double > z, complex< double > C1, complex< double > R) {
	Count += 1;
	double Modul = std::abs((z + C1 + R));
	return Modul * Modul;
}
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
double FMNGS(complex< double > x, double alpha, complex < double> Grad) {
	complex<double> GradNew;
	GradNew = Gradient(x - alpha * Grad);
	GradNew = GradNew / sqrt(real(GradNew) * real(GradNew) + imag(GradNew) * imag(GradNew));
	
	double Skalar = abs(real(GradNew) * (-real(Grad)) - imag(Grad) * imag(GradNew));
	return Skalar;
}
double FMNGS_2(complex< double > x, double alpha, complex < double> Grad, complex< double > B1, complex< double > B2, complex< double > R_1) {
	complex<double> GradNew;
	GradNew = Gradient_2(x - alpha * Grad, B1, B2, R_1);
	GradNew = GradNew / sqrt(real(GradNew) * real(GradNew) + imag(GradNew) * imag(GradNew));

	double Skalar = abs(real(GradNew) * (-real(Grad)) - imag(Grad) * imag(GradNew));
	return Skalar;
}
double FMNGS_3(complex< double > x, double alpha, complex < double> Grad, complex< double > A1, complex< double > R_2) {
	complex<double> GradNew;
	GradNew = Gradient_3(x - alpha * Grad, A1, R_2);
	GradNew = GradNew / sqrt(real(GradNew) * real(GradNew) + imag(GradNew) * imag(GradNew));

	double Skalar = abs(real(GradNew) * (-real(Grad)) - imag(Grad) * imag(GradNew));
	return Skalar;
}





double GoldenRatio(complex<double> x, double alpha, complex<double> Grad){
	/*cout << "Golden ratio method" << endl;*/
	double a = 0.00000000000000001;
	double b = alpha;
	double epsilon = 0.000001;
	
	double stopper = (b - a) / 2;
	double FC, FD;
	FC = FD = 0;
	double c, d;
	cout << setw(9) << "Iteration " << setw(12) << "a "
		<< setw(9) << "b" << setw(12) << "Stopper"
		<< setw(9) << "c" << setw(12) << "F(c)"
		<< setw(9) << "d" << setw(12) << "F(d)" << endl;
	while (stopper > epsilon) {

		c = ((3 - sqrt(5)) / 2) * (b - a) + a;
		d = ((sqrt(5) - 1) / 2) * (b - a) + a;
		if (FC == 0) {
			FC = Func(x - c * Grad, C1,C2,C3);
		}
		if (FD == 0) {
			FD = Func(x - d * Grad, C1, C2, C3);
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
		
	};
	/*cout << "Min of golden ratio = " << setprecision(4) << (a + b) / 2 << endl;*/
	return (a + b) / 2;
}
double GoldenRatio_2(complex<double> x, double alpha, complex<double> Grad, complex<double> B1, complex<double> B2, complex<double> R_1) {
	/*cout << "Golden ratio method" << endl;*/
	double a = 0.00000000000000001;
	double b = alpha;
	double epsilon = 0.000001;
	iteration = 1;
	double stopper = (b - a) / 2;
	double FC, FD;
	FC = FD = 0;
	double c, d;
	/*cout << setw(9) << "Iteration " << setw(12) << "a "
		<< setw(9) << "b" << setw(12) << "Stopper"
		<< setw(9) << "c" << setw(12) << "F(c)"
		<< setw(9) << "d" << setw(12) << "F(d)" << endl;*/
	while (stopper > epsilon) {

		c = ((3 - sqrt(5)) / 2) * (b - a) + a;
		d = ((sqrt(5) - 1) / 2) * (b - a) + a;
		if (FC == 0) {
			FC = Func_2(x - c * Grad, B1, B2, R_1);
		}
		if (FD == 0) {
			FD = Func_2(x - d * Grad, B1, B2, R_1);
		}
		/*cout << setw(10) << setprecision(6) << iteration << setw(13) << setprecision(6) << a
			<< setw(10) << setprecision(6) << b << setw(13) << setprecision(6) << stopper
			<< setw(10) << setprecision(6) << c << setw(13) << setprecision(6) << FC
			<< setw(10) << setprecision(6) << d << setw(13) << setprecision(6) << FD << endl;*/


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
	/*cout << "Min of golden ratio = " << setprecision(4) << (a + b) / 2 << endl;*/
	return (a + b) / 2;
}
double GoldenRatio_3(complex<double> x, double alpha, complex<double> Grad, complex<double> A1, complex<double> R_2) {
	/*cout << "Golden ratio method" << endl;*/
	double a = 0.00000000000000001;
	double b = alpha;
	double epsilon = 0.000001;
	iteration = 1;
	double stopper = (b - a) / 2;
	double FC, FD;
	FC = FD = 0;
	double c, d;
	/*cout << setw(9) << "Iteration " << setw(12) << "a "
		<< setw(9) << "b" << setw(12) << "Stopper"
		<< setw(9) << "c" << setw(12) << "F(c)"
		<< setw(9) << "d" << setw(12) << "F(d)" << endl;*/
	while (stopper > epsilon) {

		c = ((3 - sqrt(5)) / 2) * (b - a) + a;
		d = ((sqrt(5) - 1) / 2) * (b - a) + a;
		if (FC == 0) {
			FC = Func_3(x - c * Grad, A1, R_2);
		}
		if (FD == 0) {
			FD = Func_3(x - d * Grad, A1, R_2);
		}
		/*cout << setw(10) << setprecision(6) << iteration << setw(13) << setprecision(6) << a
			<< setw(10) << setprecision(6) << b << setw(13) << setprecision(6) << stopper
			<< setw(10) << setprecision(6) << c << setw(13) << setprecision(6) << FC
			<< setw(10) << setprecision(6) << d << setw(13) << setprecision(6) << FD << endl;*/


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
	/*cout << "Min of golden ratio = " << setprecision(4) << (a + b) / 2 << endl;*/
	return (a + b) / 2;
}
double GR(complex<double> x, double alpha, complex<double> Grad, double F)	
{
	double minimum, F_New, Fm, shag;
	minimum = 0;
	Fm = F;
	shag = alpha;
	while (shag > 0.0001) {
		while (shag > 0.0001) {
			F_New = Func(x - shag * Grad, C1, C2, C3);
			if (F_New > Fm) {
				shag *= 0.51;
				break;
			}
			else {
				Fm = F_New;
				x = x - shag * Grad;
				minimum += shag;
			}
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x - shag * Grad
				<< setw(30) << F_New << setw(30) << Grad << setw(30) << shag
				<< endl;
		};
		while (shag > 0.0001) {
			F_New = Func(x + shag * Grad, C1, C2, C3);
			if (F_New > Fm) {
				shag *= 0.51;
				break;
			}
			else {
				Fm = F_New;
				x = x - shag * Grad;
				minimum -= shag;
			}
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x + shag * Grad
				<< setw(30) << F_New << setw(30) << Grad << setw(30) << shag
				<< endl;
		}
	}
	
	return minimum;
}
double GR_2(complex<double> x, double alpha, complex<double> Grad, complex<double> B1, complex<double> B2, complex<double> R_1)
{
	double minimum, F_New, F, shag;
	minimum = 0;
	F = 1231231212132;
	shag = alpha;
	while (shag > 0.0001) {
		while (shag > 0.00001) {
			F_New = Func_2(x - shag * Grad, B1, B2, R_1);
			if (F_New > F) {
				shag *= 0.51;
				break;
			}
			else {
				F = F_New;
				x = x - shag * Grad;
				minimum += shag;
			}
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x - shag * Grad
				<< setw(30) << F_New << setw(30) << Grad << setw(30) << shag
				<< endl;
		};
		while (shag > 0.00001) {
			F_New = Func_2(x + shag * Grad, B1, B2, R_1);
			if (F_New > F) {
				shag *= 0.51;
				break;
			}
			else {
				F = F_New;
				x = x - shag * Grad;
				minimum -= shag;
			}
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x + shag * Grad
				<< setw(30) << F_New << setw(30) << Grad << setw(30) << shag
				<< endl;
		}
	}
	return minimum;
}
double GR_3(complex<double> x, double alpha, complex<double> Grad, complex<double> A1, complex<double> R_2)
{
	double minimum, F_New, F, shag;
	minimum = 0;
	F = 1231231212132;
	shag = alpha;
	while (shag > 0.0001) {
		while (shag > 0.0001) {
			F_New = Func_3(x - shag * Grad, A1, R_2);
			if (F_New > F) {
				shag *= 0.51;
				break;
			}
			else {
				F = F_New;
				x = x - shag * Grad;
				minimum += shag;
			}
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x - shag * Grad
				<< setw(30) << F_New << setw(30) << Grad << setw(30) << shag
				<< endl;
		};
		while (shag > 0.0001) {
			F_New = Func_3(x + shag * Grad, A1, R_2);
			if (F_New > F) {
				shag *= 0.51;
				break;
			}
			else {
				F = F_New;
				x = x - shag * Grad;
				minimum -= shag;
			}
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x + shag * Grad
				<< setw(30) << F_New << setw(30) << Grad << setw(30) << shag
				<< endl;
		}
	}
	return minimum;
}
int main() {
	complex<double> x(6, 9);
	complex<double> Grad(0, 0);
	complex<double> NewX(0, 1);
	complex<double> GradNew(0, 0);
	complex<double> GradOut(0, 0);
	complex<double> X_new;
	double F, alpha, cosinus, start_alpha, F_new;
	
	iteration = 0;
	GradOut = Grad = Gradient(x);
	double norm = sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	
	start_alpha = alpha = 4;
	cout << setw(8) << "Iteration " << setw(30) << "x"
		<< setw(30) << "F" << setw(30) << "Grad" << setw(30) << "alpha" << setw(30) << "COS" << setw(30) << "FMNGS"
		 << endl;
	F_new = 1231213321231;
	while (alpha > 0.0005) {

		F = Func(x, C1, C2, C3);
		GradOut = Grad;
		X_new = x - alpha * Grad;
		F_new = Func(X_new, C1, C2, C3);
		std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x
			<< setw(30) << F << setw(30) << Grad << setw(30) << alpha <<  setw(30) << "START OF MIN"
			<< endl;
		std::cout << setw(8) << iteration << setw(30) << setprecision(8) << X_new
			<< setw(30) << F_new << setw(30) << Grad << setw(30) << alpha << setw(30) << "start"
			<< endl;
		/*alpha = GR(x, alpha, Grad, F);*/
		while (F_new > F) {
			alpha *= 0.5;
			X_new = x - alpha * Grad;
			F_new = Func(X_new, C1, C2, C3);
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << X_new
				<< setw(30) << F_new << setw(30) << Grad << setw(30) << alpha 
				<< endl;
		}
		
		
		alpha *= 2;
		alpha = GoldenRatio(x, alpha, Grad);
		EPS *= 0.1;
		
		x = x - alpha * Grad;
		F = Func(x,C1,C2,C3);
		GradNew = Gradient(x);
		Count -= 2;
		GradNew = GradNew / sqrt(real(GradNew) * real(GradNew) + imag(GradNew) * imag(GradNew));
		cosinus = (real(Grad) * real(GradNew) + imag(Grad) * imag(GradNew));
		
		Count -= 2;
		
		Grad = Gradient(x);
		norm = sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
		Grad = Grad / norm;
		
		std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x
			<< setw(30) << F << setw(30) << Grad << setw(30) << alpha << setw(30) << cosinus << setw(30) << "end"
			<< endl;
		iteration += 1;
	}std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x
		<< setw(30) << F << setw(30) << Grad << setw(30) << alpha << setw(30) << cosinus << setw(30) << "end"
		<< endl;
	Count -= 2;
	cout << "First root is " << setprecision(8) << x << endl;
	cout << "Number  of function calls: " << Count << endl;
	Count = 0;
	// Horner Scheme
	complex< double > B1, B2, R_1;
	B1 = x + -C1 - C2 - C3;
	B2 = x * B1 + (C1 * C2 + C2 * C3 + C3 * C1);
	R_1 = B2 * x + (-C1 * C2 * C3);
	cout << "Function after Horner: |z^2 + z * " << B1 << " + " << B2 << "+ " << R_1 << "|^2" << endl << endl;
	complex<double> xx(6, 9);
	x = xx;
	iteration = 0;
	GradOut = Grad = Gradient_2(x, B1, B2, R_1);
	norm = sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	F = Func_2(x, B1, B2, R_1);
	start_alpha = alpha = 10;
	cout << setw(8) << "Iteration " << setw(30) << "x"
		<< setw(30) << "F" << setw(30) << "Grad" << setw(30) << "alpha" << setw(30) << "COS" << setw(30) << "FMNGS"
		<< endl;

	while (alpha > 0.00005) {

		F = Func_2(x, B1,B2,R_1);
		GradOut = Grad;
		X_new = x - alpha * Grad;
		F_new = Func_2(X_new, B1, B2, R_1);
		std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x
			<< setw(30) << F << setw(30) << Grad << setw(30) << alpha << setw(30) << "start"
			<< endl;
		/*alpha = GR(x, alpha, Grad, F);*/
		while (F_new > F) {
			alpha *= 0.9;
			X_new = x - alpha * Grad;
			F_new = Func_2(X_new, B1, B2, R_1);
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << X_new
				<< setw(30) << F_new << setw(30) << Grad << setw(30) << alpha
				<< endl;
		}
		GradNew = Gradient_2(x - alpha * Grad, B1, B2, R_1);
		GradNew = GradNew / sqrt(real(GradNew) * real(GradNew) + imag(GradNew) * imag(GradNew));
		cosinus = (real(Grad) * real(GradNew) + imag(Grad) * imag(GradNew));
		
		
		x = X_new;
		Grad = Gradient_2(x, B1, B2, R_1);
		norm = sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
		Grad = Grad / norm;
		F = F_new;
		cout << setw(8) << iteration << setw(30) << setprecision(8) << x
			<< setw(30) << F << setw(30) << Grad << setw(30) << alpha << setw(30) << cosinus << setw(30) << "end"
			<< endl;
		iteration += 1;
	}
	std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x
		<< setw(30) << F << setw(30) << Grad << setw(30) << alpha << setw(30) << cosinus << setw(30) << "end"
		<< endl;
	
	cout << "Second root is " << setprecision(8) << x << endl;
	cout << "Number  of function calls: " << Count << endl;
	Count = 0;
	// Horner Scheme
	complex< double > A1, R_2;
	A1 = x + B1;
	R_2 = x * A1 + B2 + R_1;
	cout << "Function after Horner: |z - " << A1 << "  + " << R_2 << "|^2" << endl;
	complex<double> xxx(10, 9);
	x = xxx;
	iteration = 0;
	GradOut = Grad = Gradient_3(x, A1, R_2);
	norm = sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	Grad = Grad / sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
	F = Func_3(x, A1, R_2);
	start_alpha = alpha = 10;
	cout << setw(8) << "Iteration " << setw(30) << "x"
		<< setw(30) << "F" << setw(30) << "Grad" << setw(30) << "alpha" << setw(30) << "COS" << setw(30) << "FMNGS"
		<< endl;

	while (alpha > 0.0005) {

		F = Func_3(x, A1, R_2);
		GradOut = Grad;
		X_new = x - alpha * Grad;
		F_new = Func_3(X_new, A1, R_2);
		std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x
			<< setw(30) << F << setw(30) << Grad << setw(30) << alpha << setw(30) << "start"
			<< endl;
		/*alpha = GR(x, alpha, Grad, F);*/
		while (F_new > F) {
			alpha *= 0.9;
			X_new = x - alpha * Grad;
			F_new = Func_3(X_new, A1, R_2);
			std::cout << setw(8) << iteration << setw(30) << setprecision(8) << X_new
				<< setw(30) << F_new << setw(30) << Grad << setw(30) << alpha
				<< endl;
		}
		GradNew = Gradient_3(x - alpha * Grad, A1, R_2);
		GradNew = GradNew / sqrt(real(GradNew) * real(GradNew) + imag(GradNew) * imag(GradNew));
		cosinus = (real(Grad) * real(GradNew) + imag(Grad) * imag(GradNew));
		
		Count -= 2;
		x = X_new;
		Grad = Gradient_3(x, A1, R_2);
		norm = sqrt(real(Grad) * real(Grad) + imag(Grad) * imag(Grad));
		Grad = Grad / norm;
		F = F_new;
		cout << setw(8) << iteration << setw(30) << setprecision(8) << x
			<< setw(30) << F << setw(30) << Grad << setw(30) << alpha << setw(30) << cosinus << setw(30) << "end"
			<< endl;
		iteration += 1;
	}
	std::cout << setw(8) << iteration << setw(30) << setprecision(8) << x
		<< setw(30) << F << setw(30) << Grad << setw(30) << alpha << setw(30) << cosinus << setw(30) << "end"
		<< endl;
	Count -= 2;
	cout << "Third root is " << setprecision(8) << x << endl;
	cout << "Number  of function calls: " << Count << endl;
	return 0;
}

