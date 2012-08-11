#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "filter.h"
#include "filtcoef.h"

using namespace std;

void biorfilt(string name,Laurent<double> &lpd,Laurent<double> &hpd,Laurent<double> &lpr,
                            Laurent<double> &hpr) {
	vector<double> lp1,lp2,hp1,hp2;
	filtcoef(name,lp1,hp1,lp2,hp2);
	int pow1=0;
	//reverse(lp2.begin(),lp2.end());
	reverse(hp2.begin(),hp2.end());
	
	int i=0;
	
	while  (lp2[i] == 0.0) {
		lp2.erase(lp2.begin());
		
	}
	
	i=lp2.size();
	
	while (lp2[i-1] == 0.0) {
		lp2.erase(lp2.end()-1);
		i=lp2.size();
	}
	
	i=0;
	
	while (hp2[i] == 0.0) {
		hp2.erase(hp2.begin());
		
	}
	
	i=hp2.size();
	
	while(hp2[i-1] == 0.0) {
		hp2.erase(hp2.end()-1);
		i=hp2.size();
	}
	
	lpr.setPoly(lp2,pow1);
	hpr.setPoly(hp2,pow1);
	
	Laurent<double> mn1,mn2,flv1,flv2,pr,t1,t2,flv3,flv4;
	
	//double coeff= pow(-1.0,-1);
	vector<double> coeff_mn;
	coeff_mn.push_back(-1.0);
	mn1.setPoly(coeff_mn,-1);
	

	vector<double> coeff_mn2;
	coeff_mn2.push_back(1.0);
	mn2.setPoly(coeff_mn2,-1);
	
	flv1.nzinv(hpr);
	lpd.LaurentMult(mn1,flv1);
	
	flv2.nzinv(lpr);
	hpd.LaurentMult(mn2,flv2);

	
	flv3.zinv(lpd);
	flv4.zinv(hpd);
	
	t1.LaurentMult(flv3,lpr);
	t2.LaurentMult(flv4,hpr);
	pr.LaurentAdd(t1,t2);
		
	int inc = 0;
	
	if (pr.isMono()) {
		inc = -1 * pr.monoDeg();
	} else {
		cout << "Filters Don't Specify PR Condition" << endl;
	}
	
	string str2="bior3.";
	
	if (name.find(str2,0) != string::npos) {
		inc++;
	}
		
	hpr.setPoly(hp2,inc);
	
	
	flv1.nzinv(hpr);
	lpd.LaurentMult(mn1,flv1);

	
}

void orthfilt(string name,Laurent<double> &lpd,Laurent<double> &hpd,Laurent<double> &lpr,
                            Laurent<double> &hpr) {
	vector<double> lp1,lp2,hp1,hp2;
	filtcoef(name,lp1,hp1,lp2,hp2);
	//int N=lp1.size()-1;
	reverse(lp1.begin(),lp1.end());
	int pow1=0;
	
	lpd.setPoly(lp1,pow1);
	lpr.setPoly(lp1,pow1);
	
	Laurent<double> flp,mn;
	flp.onzinv(lpd);
	double coeff= -1.0;
	vector<double> coeff_mn;
	coeff_mn.push_back(coeff);
	mn.setPoly(coeff_mn,-1);
	
	hpd.LaurentMult(mn,flp);
	hpr=hpd;
	
}

void lpoly(string name,Laurent<double> &lpd,Laurent<double> &hpd,Laurent<double> &lpr,
                            Laurent<double> &hpr) {
	
	string fname;
    fname=name.substr(0,2);
    
    if (fname == "db" || fname == "sy" || fname == "co") {
		orthfilt(name,lpd,hpd,lpr,hpr);
	} else if (fname == "bi") {
		biorfilt(name,lpd,hpd,lpr,hpr);
	}					
								
	
	}	