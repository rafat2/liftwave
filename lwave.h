#ifndef LWAVE_H
#define LWAVE_H
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "lift.h"
#include "alg.h"

using namespace std;

template <typename T>
class lwt {
	vector<T> cA,cD;
	vector<int> cD_length;
	int level;

public:
    lwt(vector<T> &signal, liftscheme &lft){
	level=1;	
	vector<double> coeff;
	vector<int> lenv;
	string lat;
	double K;
	lft.getScheme(coeff,lenv,lat,K);
	
	// Number Of Liftin Stages N
	int N = lat.size();
	vector<T> sl,dl;
	split(signal,sl,dl);
	int cume_coeff=0;
	
	for (int i=0; i < N ; i++) {
		char lft_type = lat.at(i);
		vector<double> filt;
		int len_filt = lenv[2*i];
		int max_pow = lenv[2*i+1];
		
		for (int j=0; j < len_filt; j++) {
			filt.push_back(coeff[cume_coeff+j]);
		}
		cume_coeff=cume_coeff+len_filt;
		
		if (lft_type == 'd') {
			
			for (int len_dl = 0; len_dl < (int) dl.size();len_dl++) {
				T temp = (T) 0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_dl+max_pow-lf) >= 0 && (len_dl+max_pow-lf) < (int) sl.size()) {
						temp=temp+filt[lf]*sl[len_dl+max_pow-lf];
					
					}
				}
				dl[len_dl]=dl[len_dl]-temp;
				
			}
			
		} else if (lft_type == 'p') {
			
			for (int len_sl = 0; len_sl < (int) sl.size();len_sl++) {
				T temp = (T) 0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_sl+max_pow-lf) >= 0 && (len_sl+max_pow-lf) < (int) dl.size()) {
						temp=temp+filt[lf]*dl[len_sl+max_pow-lf];
					
					}
				}
				sl[len_sl]=sl[len_sl]-temp;
				
			}
			
		}
		
	}
	double K1 = 1.0/K;
	vecmult(sl,(T) K1);
	vecmult(dl,(T) K);
	cA=sl;
	cD=dl;
	cD_length.clear();
	cD_length.push_back((int) cD.size());
		
	}
	
	lwt(vector<T> &signal, string &name){
	level=1;	
	liftscheme lft(name);
	lwt<T> wavelift(signal,lft);
	vector<T> sx,dx;
	wavelift.getCoeff(sx,dx);
	cA=sx;
	cD=dx;
	cD_length.clear();
	cD_length.push_back((int) cD.size());
	}
	
	lwt(vector<T> &signal, liftscheme &lft, int &J) {
	/*	int Max_Iter;
		Max_Iter = (int) ceil(log( double(signal.size()))/log (2.0)) - 1;

		if ( Max_Iter < J) {
			J = Max_Iter;

		}*/
		
		vector<T> temp=signal;
		vector<T> det,temp2;
		vector<int> len_det;
		
		for (int iter=0; iter < J; iter++) {
			lwt jlevel(temp,lft);
			jlevel.getCoeff(temp,temp2);
			int len_d = temp2.size();
			det.insert(det.begin(),temp2.begin(),temp2.end());
			len_det.insert(len_det.begin(),len_d);
		}
		cA=temp;
		cD=det;
		cD_length=len_det;
		level = J;
		
	}
	
	lwt(vector<T> &signal, string &name, int &J){
	liftscheme lft(name);
	lwt<T> wavelift(signal,lft);
	vector<T> sx,dx;
	wavelift.getCoeff(sx,dx);
	cA=sx;
	cD=dx;
	vector<int> cdlen;
	wavelift.getDetailVec(cdlen);
	cD_length=cdlen;
	level=J;
	}
	
void getCoeff(vector<T> &appx, vector<T> &det) {
	appx = cA;
	det = cD;
}

void getDetailVec(vector<int> &detvec) {
	detvec=cD_length;
}

int getLevels() {
	return level;
}
	
	
	virtual ~lwt()
	{
	}

};

template <typename T>
class ilwt {
	vector<T> signal;

public:
	ilwt(vector<T> &sl, vector<T> &dl, liftscheme &lft){
	vector<double> coeff;
	vector<int> lenv;
	string lat;
	double K;
	lft.getScheme(coeff,lenv,lat,K);
	//vector<T> sl,dl;
	//sl=cA;
	//dl=cD;
	double K1 = 1.0/K;
	vecmult(sl,(T) K);
	vecmult(dl,(T) K1);
	
	int N = lat.size();
	
	int cume_coeff=coeff.size();
	
	for (int i=N-1; i >= 0 ; i--) {
		char lft_type = lat.at(i);
		vector<double> filt;
		int len_filt = lenv[2*i];
		int max_pow = lenv[2*i+1];
		
		cume_coeff=cume_coeff-len_filt;
		
		for (int j=0; j < len_filt; j++) {
			filt.push_back(coeff[cume_coeff+j]);
		}
		
		
		if (lft_type == 'd') {
			
			for (int len_dl = 0; len_dl < (int) dl.size();len_dl++) {
				T temp = (T) 0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_dl+max_pow-lf) >= 0 && (len_dl+max_pow-lf) < (int) sl.size()) {
						temp=temp+filt[lf]*sl[len_dl+max_pow-lf];
					
					}
				}
				dl[len_dl]=dl[len_dl]+temp;
				
			}
			
		} else if (lft_type == 'p') {
			
			for (int len_sl = 0; len_sl < (int) sl.size();len_sl++) {
				T temp = (T) 0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_sl+max_pow-lf) >= 0 && (len_sl+max_pow-lf) < (int) dl.size()) {
						temp=temp+filt[lf]*dl[len_sl+max_pow-lf];
					
					}
				}
				sl[len_sl]=sl[len_sl]+temp;
				
			}
			
		}
		
	}
	vector<T> idwt_oup;
	merge(idwt_oup,sl,dl);
	
	signal=idwt_oup;
	
	}
	
	ilwt(vector<T> &sl,vector<T> &dl, string &name){
	liftscheme lft(name);
	ilwt<T> wavelift(sl,dl,lft);
	vector<T> sigx;
	wavelift.getSignal(sigx);
	signal=sigx;
	}
	
	ilwt(lwt<T> &wt,liftscheme &lft) {
	int J=wt.getLevels();
	vector<T> sl,dl;
	wt.getCoeff(sl,dl);
	vector<int> detv;
	
	wt.getDetailVec(detv);
	int total=0;
	
	for (int i=0; i < J; i++) {
		vector<T> temp,temp2;
		for (int j=0; j < (int) detv[i]; j++) {
			temp.push_back(dl[total+j]);
		}
		total=total+(int) detv[i];
		ilwt<T> iwt(sl,temp,lft);
		iwt.getSignal(temp2);
		sl=temp2;
		
	}
	signal=sl;
	}
	
	ilwt(lwt<T> &wt, string &name){
	liftscheme lft(name);
	ilwt<T> wavelift(wt,lft);
	vector<T> sigx;
	wavelift.getSignal(sigx);
	signal=sigx;
	}
	
void getSignal(vector<T> &sig) {
	sig=signal;
}
	
	virtual ~ilwt()
	{
	}
};

template <typename T>
class lwt2 {
	vector<T> cLL,cLH,cHL,cHH;
	int rowLL,colLL;
	int rowLH,colLH;
	int rowHL,colHL;
	int rowHH,colHH;
	
public:
    lwt2(vector<T> &signal,int rows, int cols, liftscheme &lft) {
		vector<T> L,H;
		int rows_L,cols_L,rows_H,cols_H;
		rows_L=rows;
		rows_H=rows;
		for (int i=0; i < rows; i++) {
			vector<T> temp;
			temp.assign(signal.begin()+i*cols,signal.begin()+(i+1)*cols);
			lwt<T> lwt1(temp,lft);
			vector<T> a,d;
	        lwt1.getCoeff(a,d);
			L.insert(L.end(),a.begin(),a.end());
			H.insert(H.end(),d.begin(),d.end());
			if (i==0) {
				cols_L=a.size();
				cols_H=d.size();
			}
			
		}
		
		vector<T> LT,HT;
		transpose(L,rows_L,cols_L,LT);
		transpose(H,rows_H,cols_H,HT);
		int rows_ll,cols_ll,rows_lh,cols_lh;
		
		vector<T> LL,LH;
		
		// Low Pass Stage
		cols_ll=cols_L;
		cols_lh=cols_L;

		for (int i=0; i < cols_L; i++) {
			vector<T> temp;
			temp.assign(LT.begin()+i*rows_L,LT.begin()+(i+1)*rows_L);
			lwt<T> lwt1(temp,lft);
			vector<T> a,d;
	        lwt1.getCoeff(a,d);
			LL.insert(LL.end(),a.begin(),a.end());
			LH.insert(LH.end(),d.begin(),d.end());
			if (i==0) {
				rows_ll=a.size();
				rows_lh=d.size();
			}
			
		}
		
		
		int rows_hl,cols_hl,rows_hh,cols_hh;
		
		vector<T> HL,HH;
		
		// High Pass Stage
		cols_hl=cols_H;
		cols_hh=cols_H;
		for (int i=0; i < cols_H; i++) {
			vector<T> temp;
			temp.assign(HT.begin()+i*rows_H,HT.begin()+(i+1)*rows_H);
			lwt<T> lwt1(temp,lft);
			vector<T> a,d;
	        lwt1.getCoeff(a,d);
			HL.insert(HL.end(),a.begin(),a.end());
			HH.insert(HH.end(),d.begin(),d.end());
			if (i==0) {
				rows_hl=a.size();
				rows_hh=d.size();
			}
			
		}
		
		
		
		//cLL=LL;cLH=LH;cHL=HL;cHH=HH;
		transpose(LL,cols_ll,rows_ll,cLL);
		transpose(LH,cols_lh,rows_lh,cLH);
		transpose(HL,cols_hl,rows_hl,cHL);
		transpose(HH,cols_hh,rows_hh,cHH);
		rowLL=rows_ll;rowLH=rows_lh;
		rowHL=rows_hl;rowHH=rows_hh;
		colLL=cols_ll;colLH=cols_lh;
		colHL=cols_hl;colHH=cols_hh;
		
		
	}	
	
void getCoef(vector<T> &aLL, vector<T> &aLH, vector<T> &aHL, vector<T> &aHH) {
	aLL=cLL;
	aLH=cLH;
	aHL=cHL;
	aHH=cHH;
	
}	

void getDim(vector<int> &dimvec) {
	dimvec.push_back(rowLL);
	dimvec.push_back(colLL);
	dimvec.push_back(rowLH);
	dimvec.push_back(colLH);
	dimvec.push_back(rowHL);
	dimvec.push_back(colHL);
	dimvec.push_back(rowHH);
	dimvec.push_back(colHH);
}
	
	virtual ~lwt2() {
		
	}
};

template <typename T>
class ilwt2 {
	
public:
    ilwt2(lwt2<T> &wt,liftscheme &lft) {
		
		
	}	
	
	virtual ~ilwt2() {
		
	}
};

#endif // LWAVE_H
