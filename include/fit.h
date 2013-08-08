#ifndef FIT_H_
#define FIT_H_

#include <complex>
#include <string>
typedef std::complex<double> COMPLEX; 


static const double ascale = 0.1209;
static const double hbarc = 197.327;
static const double ToMev = hbarc/ascale;

//set in out info

//#define m_po 0.597  // redused mass of PROTON_omega Mev
#define binPot(ixyz,iconf) binPot[ixyz + XYZnodeSites*iconf]



void fit(double* datain,double* datain_sigma, double& b1, double& b2,double& b3,double& b4,double& b5,double& b6,double& b7,double& b8,double& b9,double& chisq);
void fit(double* datain,double* datain_sigma, double& b1, double& b2,double& b3,double& b4,double& b5,double& b6,double& chisq);
void fit(double* datain,double* datain_sigma, double& b1, double& b2,double& b3,double& chisq);

inline void reScale(COMPLEX* pot){
  for(int id = 0; id < XYZnodeSites; id++){
    pot[id] = ToMev*pot[id];
  }
}

extern "C" {
  void mrqmin_(double x[],double y[],double sig[], int& ndata,double a[],int ia[],int& ma,
	       double covar[], double alpha[],int& nca, double &chisq,
	       void (*func)(double& x,double a[],double& yfit,double dyda[],int& ma),double& alambda);	       
  void single_exp(double& x,double a[],double& yfit,double dyda[],int& ma);
}



#endif
