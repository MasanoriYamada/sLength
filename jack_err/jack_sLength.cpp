//--------------------------------------------------------------------------
/**
 * @Filex main_fit_pot3g.cpp
 * @brief fiting potential 3 gauss
 * @ingroup YAMADA
 * @author  M.YAMADA * @date    Sat Jun 13 22:09:45 2013
 */
//--------------------------------------------------------------------------


#include <new>
#include "../include/io.h"
#include "../include/analys.h"
#include "../include/jack.h"
#include "../include/fit.h"


using namespace std;

std::string physInfo = "RC16x32_B1830Kud013760Ks013710C1761";
std::string inStaticsInfo = "1g1y3D_sLength";
std::string outStaticsInfo1 = "jack_1g1y3D_sL";
std::string outStaticsInfo2 = "jack_1g1y3D_effR";
std::string outPath = "/home/sinyamada/results/set1/Spin0-0Bin50/sLength/jack_error";




bool inBinary = false;
bool outBinary = false; 
const int dataSize = 3;
const int E_in = 1;
const int E_Max = 30;
const int NE = E_Max +1;

main(){


  IODATA Data;

  Data.setReadBinaryMode(inBinary);
  Data.setWriteBinaryMode(outBinary);
  Data.setConfSize(binnumber);
  Data.setD("1D");

  root_mkdir(outPath.c_str());
  for(int iT=T_in  ; iT< T_fi +1  ; iT++ ){
    double * ave_s = new double[NE]();
    double * err_s = new double[NE]();
    double * ave_e = new double[NE]();
    double * err_e = new double[NE]();
    double * xdata = new double[NE]();
  for (int E = E_in ;E <= E_Max ; E = E ++){
    xdata[E] = E;
    std::string inPath = "/home/sinyamada/results/set1/Spin0-0Bin50/sLength/MaxE";
    stringstream oss;
    oss<<E;
    inPath = inPath + oss.str();
    //    cout <<inPath<<endl;
    JACK *jacksLength =new JACK[NE];
    jacksLength[E].set(Confsize,binsize,dataSize);
    //    JACK jacksLength(Confsize,binsize,dataSize);

     for (int iconf=0; iconf< binnumber; iconf++) {
       double * y = new double[dataSize]();
       double * yp = new double[dataSize]();
       Data.callData(y,2,inPath,inStaticsInfo,physInfo,iconf,iT);
                     for(int id = 0;id<dataSize;id++){cout<<y[0]<<endl;}
       yp[0] = y[0];
       yp[1] = y[1];
       yp[2] = y[2];
       
       jacksLength[E].setBinData(yp,iconf);
       delete [] y;
     }
    
  double* err= new double[dataSize]();
  double* ave= new double[dataSize]();
  err= jacksLength[E].calcErr();
  ave= jacksLength[E].calcAve();
  //  for(int id = 0;id<dataSize;id++){cout<<ave[id]<<" "<<err[id]<<endl;}
  
  ave_s[E] = ave [0];
  err_s[E] = err [0];
  ave_e[E] = ave [1];
  err_e[E] = err [1];

  }//E
  //  for(int E = E_in;E<E_Max;E++){cout<<xdata[E]<<" "<<ave_s[E]<<" "<<err_s[E] <<" "<<ave_e[E]<<" "<<err_e[E]<<endl;}
  Data.outErr(xdata,ave_s,err_s,outPath,outStaticsInfo1,physInfo,0,iT,E_Max);   
  Data.outErr(xdata,ave_e,err_e,outPath,outStaticsInfo2,physInfo,0,iT,E_Max);   
    delete [] ave_s;
    delete [] err_s;
    delete [] ave_e;
    delete [] err_e;
    delete [] xdata;
  }//It
  cout <<"@End all jobs"<<endl; 
}


