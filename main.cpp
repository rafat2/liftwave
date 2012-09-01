#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "lwave.h"


using namespace std;

int main()
{
	
	string name="db2";
	int J=1;
	/*
	double lp1_a[] = {0.9501,0.2311,0.6068,0.4860,0.8913,0.7621,0.4565,0.0185,0.8214,
	0.4447,0.6154,0.7919,0.9218,0.7382,0.1763,0.4057,0.254};
	double lp1_a[] = {1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
	1.000,1.000,1.000,1.000,1.000,1.000,1.000};*/
	vector<double> sig;
    //sig.assign(lp1_a,lp1_a + sizeof(lp1_a)/sizeof(double));
	
	for (int i=0; i < 88;i ++) {
		sig.push_back(1.000);
	}
	
	int rows=11;
	int cols=8;
	
	liftscheme blift(name);
	lwt2<double> lift2(sig,rows,cols,blift);
	
	vector<double> A,B,C,D;
	lift2.getCoef(A,B,C,D);
	vector<int> lengths;
	lift2.getDim(lengths);
	cout << A.size() << " " << B.size() << " " << C.size() << " " << D.size() << endl;
	cout << lengths.size() << endl;
	for (int i=0; i < (int) B.size(); i++) {
		cout << i << " " << B[i] << endl;
	}
	
	
	 

	return 0;
}