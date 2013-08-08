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
std::string inPath = "../out";
std::string outPath = ".";
std::string physInfo = "RC16x32_B1830Kud013760Ks013710C17610.kap_013710";
std::string inStaticsInfo = "param";
std::string outStaticsInfo = "jack";
bool inBinary = false;
bool outBinary = false; 


main(){


  IODATA inPot;
  inPot.setReadBinaryMode(inBinary);
  inPot.setWriteBinaryMode(outBinary);
  inPot.setConfSize(binnumber);


  for(int iT=T_in  ; iT< T_fi +1  ; iT++ ){
    JACK jackPot(Confsize,binsize,XYZnodeSites);
    double * x = new double[XYZnodeSites]();
    for (int iconf=0; iconf< binnumber; iconf++) {
      double* a = new double[6]();
      double * y = new double[XYZnodeSites]();
      inPot.callData(a,inPath,inStaticsInfo,physInfo,iconf,iT);
#define radius_sq(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z))
#define min(a,b) (((a) < (b)) ? (a) : (b))
      //      cout <<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<" "<<a[4]<<" "<<a[5]<<endl;
      for(int iz = 0; iz<ZnodeSites; iz++){
	for(int iy = 0; iy<YnodeSites; iy++){
	  for(int ix = 0; ix<XnodeSites; ix++){
	    int ixyz = ix+ XnodeSites*(iy+YnodeSites*iz); 
	    x[ixyz] = ascale*sqrt(radius_sq( min(ix,XnodeSites - ix), min(iy,YnodeSites - iy), min(iz,ZnodeSites - iz) ));
	    y[ixyz] =a[0]*exp(-(a[1])*x[ixyz]*x[ixyz])+ a[2]*exp(-(a[3])*x[ixyz]*x[ixyz])+ a[4]*exp(-(a[5])*x[ixyz]*x[ixyz]);
	    //	    cout << x[ixyz]<<" "<<y[ixyz]<<endl;
	  }
	}
      }
      delete []a;
      delete[] y;
      jackPot.setBinData(y,iconf);
     }
    
    double* err= new double[XYZnodeSites]();
    double* ave= new double[XYZnodeSites]();
    err= jackPot.calcErr();
    ave= jackPot.calcAve();
        for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++){cout<<x[ixyz]<<" "<<ave[ixyz]<<" "<<err[ixyz]<<endl;}
    //inPot.outErr(ave,err,outPath,outStaticsInfo,physInfo,700,iT,XYZnodeSites);
        delete [] ave;
        delete [] err;
    
  }//It
  cout <<"@End all jobs"<<endl; 
}


