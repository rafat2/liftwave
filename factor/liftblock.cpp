#ifndef LIFT_H
#define LIFT_H
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "Laurent.h"

template <typename T>
class liftblock {
	int stages;
	string ltype;
	vector<Laurent<T> > lcoef;
	vector<T> Kconst;

public:
	liftblock(string &t,vector<Laurent<T> > &lc,vector<T> Kc ) {
		stages = t.size();
		ltype = t;
		lcoef = lc;
		Kconst = Kc;
	}
	
	void disp() {
		cout << "Total Number of Stages : " << stages << endl;
		cout << "--------------------------" << endl;
		for (int i=0; i < stages; i++) {
			cout << "Stage : " << i+1 << endl;
			if (ltype.at(i) == 'p') {
				cout << "Primal Lifting" << endl;
			} else if (ltype.at(i) == 'd') {
				cout << "Dual Lift" << endl;
			}
			vector<T> coeff_lc;
			lcoef[i].getPoly(coeff_lc);
            cout << "Coefficients At Stage " << i+1 << endl; 
			cout << "[ "; 
			for (int j=0; j < (int) coeff_lc.size(); j++) {
				cout << coeff_lc[j] << " ";
				
			}
			cout << " ]" << endl;
			cout << "Maximum Degree : " << lcoef[i].highdeg() << endl;
		}
	}
	
void addlift(string &c, Laurent<T> unit) {
		ltype=ltype+c;
		stages = ltype.size();
		lcoef.push_back(unit);
	}
	
virtual ~liftblock() {
		
	}

};

template <typename T>
void liftwave(string &name) {
	string pdstr;
	vector<Laurent<double> > liftcont;
	vector<double> K1K;
	
	switch(name) {
		case "db2":
		    pdstr="pdp";
		    break;
		default:
            break;	
		
		
	}
	liftblock<double> setlift(pdstr,liftcont,K1K);
	
}

#endif // LIFT_H
