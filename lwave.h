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

public:
    lwt(vector<T> &signal, liftscheme &lft){
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
					if ((len_dl+max_pow-lf) >= 0 && (len_dl+max_pow-lf) < sl.size()) {
						temp=temp+filt[lf]*sl[len_dl+max_pow-lf];
					
					}
				}
				dl[len_dl]=dl[len_dl]-temp;
				
			}
			
		} else if (lft_type == 'p') {
			
			for (int len_sl = 0; len_sl < (int) sl.size();len_sl++) {
				T temp = (T) 0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_sl+max_pow-lf) >= 0 && (len_sl+max_pow-lf) < dl.size()) {
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
		
	}
	
	lwt(vector<T> &signal, string &name){
	liftscheme lft(name);
	lwt<T> wavelift(signal,lft);
	vector<T> sx,dx;
	wavelift.getCoeff(sx,dx);
	cA=sx;
	cD=dx;
	}
	
void getCoeff(vector<T> &appx, vector<T> &det) {
	appx = cA;
	det = cD;
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
					if ((len_dl+max_pow-lf) >= 0 && (len_dl+max_pow-lf) < sl.size()) {
						temp=temp+filt[lf]*sl[len_dl+max_pow-lf];
					
					}
				}
				dl[len_dl]=dl[len_dl]+temp;
				
			}
			
		} else if (lft_type == 'p') {
			
			for (int len_sl = 0; len_sl < (int) sl.size();len_sl++) {
				T temp = (T) 0.0;
				for (int lf=0; lf < len_filt; lf++) {
					if ((len_sl+max_pow-lf) >= 0 && (len_sl+max_pow-lf) < dl.size()) {
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
	
void getSignal(vector<T> &sig) {
	sig=signal;
}
	
	virtual ~ilwt()
	{
	}
};

#endif // LWAVE_H
