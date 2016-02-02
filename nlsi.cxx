#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);

void step_lin(cmplx* const psi1, cmplx* const psi0,
          const double dt, const double dx, const int N);

void step_nichtlin(cmplx* const psi1, cmplx* const psi0,
          const double dt, const double dx, const int N);
//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;
	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-1; j++) {
			step_lin(psi1, psi0, dt/2.0, dx, Nx);
			step_nichtlin(psi0, psi1, dt, dx, Nx);
			step_lin(psi1, psi0, dt/2.0, dx, Nx);
 			h = psi0;
 			psi0 = psi1;
 			psi1 = h;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}

	return 0;
}
//-----------------------------------
void step_lin(cmplx* const psi1, cmplx* const psi0,
          const double dt, const double dx, const int N)
{
  const cmplx alpha = cmplx(0,-dt/(dx*dx));

  cmplx* d=new cmplx[N];

  for(int i=0;i<N;i++)
	{
		d[i] = 1.0 + 2.0*alpha;
	}
	
  for(int i=1;i<N;i++)
	{
    d[i]  -=   alpha*alpha/d[i-1];
    psi0[i] -=    -alpha/d[i-1]*psi0[i-1];
  }

  psi1[N-1] = psi0[N-1]/d[N-1];
  for(int i=N-2;i>=0; i--)
    psi1[i] = (psi0[i] + alpha*psi1[i+1])/d[i];
  
	delete[] d;
}
//-----------------------------------
void step_nichtlin(cmplx* const psi1, cmplx* const psi0,
          const double dt, const double dx, const int N)
{
	for(int i=0 ; i<N ; i++)
		psi1[i] = psi0[i]*exp(cmplx(0,- norm(psi0[i])*dt));
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
//		psi0[i] = f*exp(-pow(x/10.0,2));
	}
}
