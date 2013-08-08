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


main(){


  IODATA inPot;
  inPot.setReadBinaryMode(inBinary);
  inPot.setWriteBinaryMode(outBinary);
  inPot.setConfSize(binnumber);


  for(int iT=T_in  ; iT< T_fi +1  ; iT++ ){
    JACK jackPot(Confsize,binsize,XYZnodeSites);
    COMPLEX* binPot = new COMPLEX[XYZnodeSites*binnumber]();
      double * x = new double[XYZnodeSites]();
    for (int iconf=0; iconf< binnumber; iconf++) {
      COMPLEX* pot = new COMPLEX[XYZnodeSites]();
      inPot.callData(pot,inPath,inStaticsInfo,physInfo,iconf,iT);
      reScale(pot);
      memcpy(binPot + iconf*XYZnodeSites, pot, sizeof(pot)*XYZnodeSites*2);    
      double * y = new double[XYZnodeSites]();
#define radius_sq(x,y,z) ((x)*(x) + (y)*(y) + (z)*(z))
#define min(a,b) (((a) < (b)) ? (a) : (b))

      for(int iz = 0; iz<ZnodeSites; iz++){
	for(int iy = 0; iy<YnodeSites; iy++){
	  for(int ix = 0; ix<XnodeSites; ix++){
	    int ixyz = ix+ XnodeSites*(iy+YnodeSites*iz); 
	    x[ixyz] = ascale*sqrt(radius_sq( min(ix,XnodeSites - ix), min(iy,YnodeSites - iy), min(iz,ZnodeSites - iz) ));
	    y[ixyz] = pot[(ix) +XnodeSites*((iy) + YnodeSites*((iz)))].real();
	    if(x[ixyz]==4.0*ascale || x[ixyz] ==8.0*ascale)	    cout <<x[ixyz]<<" "<<pot[ixyz].real()<<endl; 
	  }
	}
      }
      jackPot.setBinData(y,iconf);
      delete [] pot;
    }
    
    double* err= new double[XYZnodeSites]();
    double* ave= new double[XYZnodeSites]();
    err= jackPot.calcErr();
    ave= jackPot.calcAve();
    //    for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++){cout<<x[ixyz]<<" "<<ave[ixyz]<<" "<<err[ixyz]<<endl;}
    inPot.outErr(ave,err,outPath,outStaticsInfo,physInfo,700,iT,XYZnodeSites);
    
    cout << "@Start fit"<<endl;
    for (int iconf=0; iconf< binnumber; iconf++) {
      double b1,b2,b3,b4,b5,b6,b7,b8,b9,chisq;
      double* binPotIn = new double[XYZnodeSites]();
      for(int ixyz = 0;ixyz<XYZnodeSites;ixyz++){
	binPotIn[ixyz] = binPot[ixyz + iconf*XYZnodeSites].real();
      }
      //      fit(binPotIn,err,b1,b2,b3,b4,b5,b6,b7,b8,b9,chisq);
      //      cout<<iconf<<"@param "<<b1<<" "<<b1<<" "<<b2<<" "<<b3<<" "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<b8<<" "<<b9<<" chisq = "<<chisq<<endl;
      double param[9]={b1,b2,b3,b4,b5,b6,b7,b8,b9};
      delete [] binPotIn;
      inPot.outData(param,outPath,outStaticsInfo,physInfo,iconf,iT,9);
    }//iconf  
    delete [] binPot;
    //    delete [] ave;
    //    delete [] err;
    
  }//It
  cout <<"@End all jobs"<<endl; 
}


