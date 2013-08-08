//--------------------------------------------------------------------------
/**
 * @Filex main_fit_pot3g.cpp
 * @brief fiting potential 3 gauss
 * @ingroup YAMADA
 * @author  M.YAMADA * @date    Sat Jun 13 22:09:45 2013
 */
//--------------------------------------------------------------------------

#include "../include/io.h"
#include "../include/analys.h"
#include "../include/jack.h"
#include "../include/fit.h"


using namespace std;

std::string inPath = "/Users/SINYAMADA/lab/data/results/bin1/Potential/binPotential/xyz";
std::string outPath = "../out";
std::string physInfo = "RC16x32_B1830Kud013760Ks013710C17610.kap_013710";
std::string inStaticsInfo = "Potential";
std::string outStaticsInfo = "param";
bool inBinary = true;
bool outBinary = false; 

main(){


  IODATA inPot;
  inPot.setReadBinaryMode(inBinary);
  inPot.setWriteBinaryMode(outBinary);
  inPot.setConfSize(binnumber);


  for(int iT=T_in  ; iT< T_fi +1  ; iT++ ){
    JACK jackPot(Confsize,binsize,XYZnodeSites);
    COMPLEX* binPot = new COMPLEX[XYZnodeSites*binnumber];
    for (int iconf=0; iconf< binnumber; iconf++) {
      COMPLEX* pot = new COMPLEX[XYZnodeSites];
      inPot.callData(pot,inPath,inStaticsInfo,physInfo,iconf,iT);
            reScale(pot);
	    memcpy(binPot + iconf*XYZnodeSites, pot, sizeof(pot)*XYZnodeSites*2);    
      //        if(iconf == 0){      for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++){ cout<<iconf<<" "<<ixyz<<" "<<pot[ixyz]<<endl;}}
      jackPot.setBinData(pot,iconf);
      delete [] pot;
    }
    
    double* err= new double[XYZnodeSites]; 
    err= jackPot.calcErr();
    for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++){cout<<"err "<<ixyz<<" "<<err[ixyz]<<endl;}
    //  inPot.outData(err,outPath,outStaticsInfo,physInfo,700,iT,XYZnodeSites);
    
    cout << "@Start fit"<<endl;
    for (int iconf=0; iconf< binnumber; iconf++) {
      double b1,b2,b3,b4,b5,b6,b7,b8,chisq;
      double* binPotIn = new double[XYZnodeSites];
      for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++){
	binPotIn[ixyz] = binPot(ixyz,iconf).real();
      }
      fit(binPotIn,err,b1,b2,b3,b4,b5,b6,chisq);
      double param[6]={b1,b2,b3,b4,b5,b6};

      inPot.outData(param,outPath,outStaticsInfo,physInfo,iconf,iT,6);
      cout<<b1<<" "<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<chisq<<endl;
    }//iT  
 
    cout <<"@End all jobs"<<endl; 
  }

}


//---------------------fit function---------------------------------------//
  void single_exp(double& x,double a[],double& yfit,double dyda[],int& ma)
  {
    
    double b1=a[0];
    double b2=a[1];
    double b3=a[2];
    double b4=a[3];
    double b5=a[4];
    double b6=a[5];
    double b7=a[6];
    double b8=a[7];
     yfit = b1*exp(-b2*x*x)+ b3*exp(-b4*x*x)+ b5*exp(-b6*x*x);
    //   yfit = b1*exp(-b2*x*x)+ b3*exp(-b4*x*x)+ b5*exp(-b6*x*x) + b1*exp(-b2*(r_sq-(x*x)))+ b3*exp(-b4*(r_sq-(x*x)))+ b5*exp(-b6*(r_sq-(x*x))) +b1*exp(-b2*(r_sq-(x*x)))+ b3*exp(-b4*(r_sq-(x*x)))+ b5*exp(-b6*(r_sq-(x*x)));
    
    
    dyda[0] = exp(-b2*x*x);
    dyda[1] = -b1*x*x*exp(-b2*x*x);
    dyda[2] = exp(-b4*x*x);
    dyda[3] = -b3*x*x*exp(-b4*x*x);
    dyda[4] = exp(-b6*x*x);
    dyda[5] = -b5*x*x*exp(-b6*x*x);
    dyda[6] = exp(-b8*x*x);
    dyda[7] = -b7*x*x*exp(-b8*x*x);
    
  }


  void fit(double datain[],double datain_sigma[], double& b1, double& b2,double& b3,double& b4,double& b5,double& b6,double& chisq)
  {
#define datain(ix,iy,iz)  datain[ix+XnodeSites*(iy+YnodeSites*(iz))]
#define datain_sigma(ix,iy,iz)  datain_sigma[ix+Nx*(iy+Ny*(iz))]
    double x[8*8*8], y[8*8*8], dy[8*8*8];
	
    int ixyz=0;
    for(int iz = 0; iz < ZnodeSites/2; iz++){
      for(int iy = 0; iy < YnodeSites/2; iy++){
	for(int ix = 0; ix < ZnodeSites/2; ix++){	
	  x[ixyz]  = ascale*sqrt(ix*ix + iy*iy + iz*iz);

	  y[ixyz]  = datain[(ix) +XnodeSites*((iy) + YnodeSites*((iz)))];
	  dy[ixyz] = datain_sigma[(ix) +XnodeSites*((iy) + YnodeSites*((iz)))];
	  ixyz++;
	}}}//ixyz
	
    int ndata = 8*8*8;
    //			static double a[2];
    //			static int    is_initial = 1;
    //			int   ia[]={1,1}, ma =2;
    //			double covar[4],alpha[4];
    //			int    nca = 2;
    static double a[6];
    static int    is_initial = 1;
    int   ia[]={1,1,1,1,1,1}, ma =6;
    double covar[36],alpha[36];
    int    nca = 6;
    double alambda;
	
    if(is_initial==1){
      is_initial=0;
      
      a[0] =-212;
      a[1] = 2.39;
      a[2] = -374;
      a[3] =1.2;
      a[4] = 4.9;
      
      //it07
      a[0]               = -157.287  ;
      a[1]               =  2.00842 ;
      a[2]               = 201.697;
      a[3]               = 12.5821 ;
      a[4]               = 2305.73 ;
      a[5]               = 69.2551 ;
      a[6]               =  -1217.77;
      a[7]               =  69.2551   ;
      //it08
      a[0]               = -150.925   ;
      a[1]               = 2.68075  ;
      a[2]               = 515.608 ;
      a[3]               = 5.7681;
      a[4]               = 2322.14 ;
      a[5]               = 65.7699 ;
      ///a[6]               = -1201.59 ;
      //a[7]               = 65.7296  ;
      //it09
      a[0]               = -533.934  ;
      a[1]               =  2.80389 ;
      a[2]               = 486.76;
      a[3]               = 3.91707 ;
      a[4]               = 2336.8 ;
      a[5]               = 59.5949 ;
      a[6]               = -1186.94;
      a[7]               = 58.7904  ;
        
    }
    //
    // The initial call
    //
    alambda = -1.0;
    mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, single_exp, alambda);
    //
    // iterations
    //
    for(int iter=0; iter<65536; iter++){
      mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, single_exp, alambda);
      if (alambda > 1.0e64) break;
    }
    if (alambda <= 1.0e64) {
      cerr << "convergence is not achieved\n";
      exit(1);
    }
    //
    // The final call
    //
    alambda = 0.0;
    mrqmin_(x,y,dy, ndata,a,ia, ma,	covar,alpha, nca,chisq, single_exp, alambda);
    b1 = a[0];
    b2 = a[1];
    b3 = a[2];
    b4 = a[3];
    b5 = a[4];
    b6 = a[5];
    //b7 = a[6];
    //b8 = a[7];

  }
  //////////////////////////////END OF FITTING /////////////////////////////////

