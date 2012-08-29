#ifndef LIFT_H
#define LIFT_H
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "Laurent.h"

//template <typename T>
class liftscheme{
	int stages;
	string ltype;
	vector<double> lcoeff;
	vector<int> plen;
    double Kconst;
	string wname;
	
public:
	liftscheme() {
		ltype="";
		stages=0;
		Kconst=1.000;
		vector<double> lcoeff;
		vector<int> plen;
		wname="lazy";
		
	}
	
	liftscheme(string &name) {
		wname=name;
//		vector<double> coeffs;
		if (name == "haar" || name == "db1" ) {
            ltype="dp";
			stages=2;
			Kconst=0.7071;
			
			//Stage 1,2
			double d1[]={1.0000};
			double p1[]={-0.5000};
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,1,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));
			
		} else if (name == "db2") {

		    ltype="pdp";
			stages=3;
			Kconst=1.93185;
			
			//Stage 1,2,3
			double p1[]={-1.73205};
			double d1[]={0.433013,-0.0669873};
			double p2[]={1.0000};
			
			lcoeff.insert(lcoeff.begin(),p2,p2+1);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,1,1};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "db3") {

		    ltype="dpdp";
			stages=4;
			Kconst=1.9182;
			
			//Stage 1,2,3,4
			
			double d1[]={-0.412287};
			double p1[]={0.352388,-1.56514};
			double d2[]={0.492152,0.0284591};
			double p2[]={-0.38962};
			
			lcoeff.insert(lcoeff.begin(),p2,p2+1);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,0,2,1,1,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "db4") {

		    ltype="pdpdp";
			stages=5;
			Kconst=2.61312;
			
			//Stage 1,2,3,4,5
			
			double p1[]={-3.10293};
			double d1[]={0.291953,-0.0763001};
			double p2[]={-1.66253,5.19949};
			double d2[]={0.0378927,-0.00672237};
			double p3[]={0.314106};
			
			lcoeff.insert(lcoeff.begin(),p3,p3+1);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+1);
			
			int pl[]={1,0,2,0,2,2,2,-2,1,3};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior2.2") {

		    ltype="dp";
			stages=2;
			Kconst=0.707107;
			
			//Stage 1,2
			
			double d1[]={0.5,0.5};
			double p1[]={-0.25,-0.25};
			
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,2,0};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior2.4") {

		    ltype="dp";
			stages=2;
			Kconst=0.707107;
			
			//Stage 1,2
			
			double d1[]={0.5,0.5};
			double p1[]={0.046875,-0.296875,-0.296875,0.046875};
			
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+4);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,4,1};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	}
	
	}

int nlifts() {
		return stages;
	}

double K() {
	return Kconst;
}

string getName() {
	return wname;
}	

void getScheme(vector<double> &coeff, vector<int> &lenvec, string &lattice,double &Kc) {
	Kc=Kconst;
	lattice=ltype;
	coeff=lcoeff;
	lenvec=plen;
	
}
	
void disp() {
	cout << "Total Number of Stages : " << stages << endl;
	cout << "--------------------------" << endl;
	int total=0;
	for (int i=0; i < stages; i++) {
			cout << "Stage : " << i+1 << endl;
			if (ltype.at(i) == 'p') {
				cout << "Primal Lifting" << endl;
			} else if (ltype.at(i) == 'd') {
				cout << "Dual Lift" << endl;
			}
			cout << "Coefficients : ";
			int t2=0;
			for (int j=0; j < plen[2*i];j++) {
				cout << lcoeff[total+j] << " ";
				t2++;
			}
			total=total+t2;
			cout << endl;
	}
	cout << "--------------------------" << endl;
	cout << " K : " << Kconst <<endl;
}	

virtual ~liftscheme() {
	
}	
	
};

#endif // LIFT_H
