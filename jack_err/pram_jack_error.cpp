/*
 *  jack.cpp
 *  
 *
 *  Created by yamada on 1/29/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
/*information 
 t=9についてjackerrorをとる
 10conf
 16^3*32
 input ./results/proj.omega-omega1.2.7.1/s3z3.RC16x32_B1830Kud013760Ks013710C1761-1-%05d0.kap_013710.000_000_000_000_it%02d",j+41,it
 output ,/results/proj.omega-omega1.2.7.1/jack_error/s3z3.RC16x32_B1830Kud013760Ks013710C1761-1-kap_013710.000_000_000_000_it%02d",it
 常に追記モードなので注意
 チェックリスト
 読み込みは成功している
 値を全てチェックしたので絶対に間違えてない
 

 */
#include<stdio.h>
#include<stdint.h>
#include<math.h>
#include<complex>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<memory.h>
#include<cstdlib>
#include<assert.h>

static const int Numberofdata=2;							//ファイルの数(これを統計として使う)
static const int datasize=6;									//１つのファイルにあるデータの数
static const int binsize=1;									//Numberofdataをわってあまりが出ない値を使わないといけない
static const int T_in=7;
static const int T_fi=7;
static const int binnumber=(Numberofdata/binsize);
static const int switch1=1;                          //0だとbainaryから読み込む　1だと普通の数字を読み込む
#define data(id,j) data[Numberofdata*id+j]					//id:dataの数		
#define I       std::complex<double>(0.0,1.0)
using namespace std;
typedef std::complex<double> COMPLEX;



int call_file(char[],float[],double[]);
void call_data(int,float[],double[]);
void out_file(int,char[],float[],double,double);
void out_data(int,int,float[],double,double);
void jack_avesub_calc(int,double[],double[]);
double jack_ave_calc(double[]);
double jack_err_calc(double[]);



main(){
	for (int it=T_in; it < (T_fi +1); it++) {
		float rad[datasize];
		double data[datasize*Numberofdata];
	call_data(it,&(rad[0]),&(data[0]));
	for (int id=0; id < datasize; id++) {
		double jack_ave_sub[binnumber];
		double ave=0.0;
		double err=0.0;
	jack_avesub_calc(id,&(data[0]),&(jack_ave_sub[0]));
	ave=jack_ave_calc(&(jack_ave_sub[0]));
	err=jack_err_calc(&(jack_ave_sub[0]));
	out_data(id,it,&(rad[0]),ave,err);
	}}
	return 0;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void call_data(int it,float rad[datasize],double data[Numberofdata*datasize]){
	for (int j=0; j<Numberofdata; j++) {
	    char fname[200]={0};
	    double local[datasize]={0};
	    		sprintf(fname,"../out/param.000700-%06d.RC16x32_B1830Kud013760Ks013710C17610.kap_013710.it%02d",j,it);
	    //sprintf(fname,"/media/HD-PNTU2/3_2/proj.PN_BSwave/s0z0.RC16x32_B1830Kud013760Ks013710C1761-1-%05d0.kap_013710.000_000_000_000_it%02d",j,it);
		 	call_file(&(fname[0]),&(rad[0]),&(local[0]));
		for (int id=0; id < datasize; id++) {
		   data(id,j)=local[id];
	}
	}
}
/*******************************************************************************************************/
//ファイルからデータを呼び出す	no (読み込むファイルパス,読みだしたout用の配列(double)	  //
/******************************************************************************************************/
int call_file(char fname[200],float rad[datasize],double local[datasize]){
	double *proj_wave = new double[datasize];  /* 読み込みデータ格納用   */
	//string ss;
    /*
	if (switch1==0) {std::ifstream ifs( fname , ios::in | ios::binary	);			// バイナリモードで 
		if(! ifs  )
		{
			cout << fname<<"ファイルが開きません";
			return 1;
		}
	
	for(int	id=0; id<datasize; ++id)  {
	ifs>>std::setw(7)>>rad[id]>>std::setw(15)>> proj_wave[id].real()>>std::setw(15)>>proj_wave[id].imag();
			local[id]=(double)proj_wave[id].real();
			//	cout << id <<local[id]<<endl;
		return 0;	
	}}*/
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
	if (switch1==1) {std::ifstream ifs( fname , ios::in );			/* テキストモードで */
		if(! ifs  )
		{
			cout << fname<<"ファイルが開きません";
			return 1;
		}
		//		getline( ifs,ss );
		for(int	id=0; id<datasize; ++id)  {
			ifs>> rad[id]>>proj_wave[id];
			// ifs>>rad[id]>> proj_wave[id].real()>>proj_wave[id].imag();
			local[id]=(double)proj_wave[id];
						cout <<std::setw(7)<<rad[id]<<"///"<<std::setw(15)<<local[id]<<endl;
			
		}}
	delete []proj_wave;
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************************/
//call_fileを利用して、複数のファイルから呼び出したデータを配列としてまとめる (読みだしたout用の配列(double)	この中のfanameをいじると入力ファイルを変えられる。  //
/******************************************************************************************************************************************/
void out_data(int id,int it,float rad[datasize],double ave,double err){
	double local1[datasize];
	
		char fname1[200]={0};
		sprintf(fname1,"./param_ave.it%02d",it);
		out_file(id,fname1,&(rad[0]),ave,err);
}
/*******************************************************************************************************/
//結果を書き出す	no (書きだすファイルパス,平均(データの大きさ),誤差(データの大きさ)	  //
/******************************************************************************************************/
void out_file(int id, char fname[200],float rad[datasize],double ave,double err){
	std::ofstream ofs( &(fname[0]),std::ios::out | std::ios::app);
//	std::ofstream ofs( &(fname[0]),std::ios::out | std::ios::trunc);
    cout<<scientific;
    ofs<<scientific;
	ofs<< ave<<"    "<<err<<endl;
    cout<< ave<<"   "<<err<<endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************/
//jackナイフのためにconf[j]を取り除き平均をとったomega_prop[j]をだす no(入力のデータ,出力のデータ)
/*******************************************************************************************************/
void jack_avesub_calc(int id,double data[datasize*Numberofdata],double jack_ave_sub[binnumber]){
	
		double full=0.0;
	for (int j=0; j<Numberofdata; j++) {
		full=full+data(id,j);
	}
	cout <<"full"<<full<<endl;
	for (int b=0; b<binnumber; b++) {
		double subfull[Numberofdata];
		double localfull=0;
	for (int j=b*binsize; j<(b+1)*binsize; j++) {
		localfull=localfull+data(id,j);
	}
		subfull[b]=full-localfull;
		cout << subfull[b]<<"="<<full<<"-"<<localfull<<endl;
	jack_ave_sub[b]=subfull[b]/((double)Numberofdata-(double)binsize);
		cout << "jack_avesub:b"<<b<<":"<<jack_ave_sub[b]<<endl;
	}
}

/********************************************************************************************************/
//jackナイフのためにconf[j]を取り除き平均をとったomega_prop[j]をだす no(入力のデータ,出力のデータ)
/*******************************************************************************************************/
double jack_ave_calc(double jack_ave_sub[binnumber]){
	double buffer=0.0;
	double ave=0.0;
	for (int b=0; b<binnumber; b++) {
			buffer=buffer+jack_ave_sub[b];
		}
		ave=buffer/(double)binnumber;
	cout << "ave"<<ave<<endl;
	return ave;
}

/*******************************************************************************************************/
//jackナイフのためにconf[j]を取り除き平均をとったomega_prop[j]をだす no(入力のデータ,出力のデータ)
/*******************************************************************************************************/
double jack_err_calc(double jack_ave_sub[binnumber]){
	double ave1=0;
	double ave2=0;
	double err=0;
	double double_jack_ave_sub[binnumber];
		for (int b=0; b<binnumber; b++) {
		double_jack_ave_sub[b]=jack_ave_sub[b]*jack_ave_sub[b];
			cout <<"double_jack"<<double_jack_ave_sub[b]<<endl;
		}
		ave1=jack_ave_calc(&(jack_ave_sub[0]));
		ave2=jack_ave_calc(&(double_jack_ave_sub[0]));
	cout << "ave1"<<ave1<<endl;
	cout << "ave2"<<ave2<<endl;
		err= sqrt(((double)binnumber -1.0)*((ave2)-(ave1)*(ave1)));
	cout << "err"<<err<<endl;
	return err;
}
