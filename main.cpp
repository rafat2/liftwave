#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "Laurent.h"
#include "filter.h"
#include "lift.h"

using namespace std;

void rbiorfilt() {
	
}


int main()
{
	cout << "Lifting Demo" << endl;
	string name="db2";
	Laurent<double> lpd,hpd,lpr,hpr;
	lpoly(name,lpd,hpd,lpr,hpr);
	Laurent<double> leven,lodd;
	EvenOdd(lpr,leven,lodd);
	vector<Laurent<double> > loup,Q;
	Div(leven,lodd,loup);
	leven.dispPoly();
	lodd.dispPoly();
	cout <<endl;
	Laurent<double> quot,rem;
	int md;
	
	for (int i=0; i < (int) loup.size() / 2;i++) {
		quot=loup[2*i];
		rem=loup[2*i+1];
		
		quot.dispPoly();
		rem.dispPoly();
		cout << endl;
	}
	leven = lodd;
	lodd = loup[3];
	Q.push_back(loup[2]);
	loup.clear();
	Div(leven,lodd,loup);
	leven.dispPoly();
	lodd.dispPoly();
	cout <<endl;
	
	for (int i=0; i < (int) loup.size() / 2;i++) {
		quot=loup[2*i];
		rem=loup[2*i+1];
		
		quot.dispPoly();
		rem.dispPoly();
		cout << endl;
	}
	
	leven = lodd;
	lodd = loup[1];
	Q.push_back(loup[0]);
	
	Div(leven,lodd,loup);
	leven.dispPoly();
	lodd.dispPoly();
	cout <<endl;
	
	cout << lodd.isZero() << endl;
	Laurent<double> o,z;
	o.One();
	z.Zero();
	
	LaurentMat<double> Mat1,Mat2,oup;
	Mat1.SZ(Q[0]);
	Mat2.TZ(Q[1]);
	oup.MatMult(Mat1,Mat2);
	oup.dispMat();
	 
	/*
	Laurent<double> lpd,hpd,lpr,hpr;
	lpoly(name,lpd,hpd,lpr,hpr);
	string s="pd";
	vector<Laurent<double> > lcont;
	lcont.push_back(lpr);
	lcont.push_back(hpr);
	vector<double> KC;
	KC.push_back(1.0);
	KC.push_back(1.0);
	liftblock<double> els(s,lcont,KC);
	els.disp();
	string str2="p";
	els.addlift(str2,lpd);
	els.disp();
	 */ 
//	Laurent<double> even,odd;
//	EvenOdd(lpd,even,odd);
//	cout << " Polyphase test" << endl; 
//	lpd.dispPoly();
//	even.dispPoly();
//	odd.dispPoly();
//    factor(name);
	return 0;
}
