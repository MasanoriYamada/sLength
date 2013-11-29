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

std::string inPath = "/home/sinyamada/results/set1/Spin0-0Bin50/phase/bin";
std::string physInfo = "RC16x32_B1830Kud013760Ks013710C1761";
std::string inStaticsInfo = "1g1yy1D_kcotdy";
std::string outStaticsInfo = "1g1yy1D_sLength";
int dataSize;  //correspond to fiting range
const int fileDatasize = 1500;
bool inBinary = false;
bool outBinary = false; 



main(){


  IODATA inKcotd;
  inKcotd.setReadBinaryMode(inBinary);
  inKcotd.setWriteBinaryMode(outBinary);
  inKcotd.setConfSize(binnumber);

  for (int E = 1 ; E <= 30 ;E++){
    dataSize = E*100/2;
    std::string outPath = "/home/sinyamada/results/set1/Spin0-0Bin50/sLength/MaxE";
    ostringstream oss;
    oss << E;
    outPath = outPath + oss.str();
    root_mkdir(outPath.c_str());

  for(int iT=T_in  ; iT< T_fi +1  ; iT++ ){
    JACK jackKcotd;
    jackKcotd.set(Confsize,binsize,dataSize);
    double* fitInX = new double[dataSize];

    for (int iconf=0; iconf< binnumber; iconf++) {
      double* binKcotd = new double[fileDatasize];
      double* freqK = new double[fileDatasize];

      inKcotd.callData(binKcotd,2,inPath,inStaticsInfo,physInfo,iconf,iT);
      inKcotd.callData(freqK,1,inPath,inStaticsInfo,physInfo,iconf,iT);
      double* fitInY = new double[dataSize];
      for (int id = 0 ;id <dataSize; id++){
	fitInX[id] = freqK[id];
	fitInY[id] = binKcotd[id];
      }
      if(iconf == 0){      for(int id = 0;id<dataSize;id++){ cout<<iconf<<" "<<fitInX[id]<<" "<<fitInY[id]<<endl;}}
      jackKcotd.setBinData(fitInY,iconf);
      delete [] fitInY;
      delete [] binKcotd;
      delete [] freqK;
    }
    
    double* err= new double[dataSize];
    double* ave= new double[dataSize]; 
    err= jackKcotd.calcErr();
    ave= jackKcotd.calcAve();
        cout << "@Start Ave Err"<<endl;
	//for(int id = 0;id<dataSize;id++){cout<<fitInX[id]<<" "<<ave[id]<<" "<<err[id]<<endl;}
    
    cout << "@Start fit"<<endl;
    for (int iconf=0; iconf< binnumber; iconf++) {
      double* binKcotd = new double[fileDatasize];
      inKcotd.callData(binKcotd,2,inPath,inStaticsInfo,physInfo,iconf,iT);
      double b1,b2,chisq;
      double* fitInY = new double[dataSize];
      for (int id = 0 ;id <dataSize; id++){
	fitInY[id] = binKcotd[id];
      }
      fit(fitInX,fitInY,err,b1,b2,chisq);
      delete [] fitInY;
      double param[3]={1/b1,b2,chisq};
      
      cout<<1/b1<<" "<<b2<<" "<<chisq<<endl;
      inKcotd.outData(param,outPath,outStaticsInfo,physInfo,iconf,iT,3);

      delete [] binKcotd;
      
    }//iconf
      delete [] fitInX; 
  }
  }
    cout <<"@End all jobs"<<endl; 

}


//---------------------fit function---------------------------------------//
  void single_exp(double& x,double a[],double& yfit,double dyda[],int& ma)
  {
    
    double b1=a[0];
    double b2=a[1];
    // for it = 7
    
    /*        yfit = 1.0/b1 + (1.0/2.0)*b2*x*x;
    
    dyda[0] = -1.0/(b1*b1);
    dyda[1] = (1.0/2.0)*x*x;
    */
          // for it = 8
        yfit = b1 + (1.0/2.0)*b2*x*x;
    
    dyda[0] = 1.0;
    dyda[1] = (1.0/2.0)*x*x;
    
  }


void fit(double x[],double y[],double dy[], double& b1, double& b2,double& chisq)
  {	
    int ndata = dataSize;
    //			static double a[2];
    //			static int    is_initial = 1;
    //			int   ia[]={1,1}, ma =2;
    //			double covar[4],alpha[4];
    //			int    nca = 2;
    static double a[2];
    static int    is_initial = 1;
    int   ia[]={1,1}, ma =2;
    double covar[4],alpha[4];
    int    nca = 2;
    double alambda;
	
    if(is_initial==1){
       is_initial=0;
      //it 7      
       // a[0] = -1.0/50.0;
       //a[1] = 0.1;
      //it 8,9
        a[0] = 0.002;
       a[1] = 0.006;
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

  }
  //////////////////////////////END OF FITTING /////////////////////////////////

