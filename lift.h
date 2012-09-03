//============================================================================
// Name : 1D/2D Wavelet Transform
// Author : Rafat Hussain
// Version :
// Copyright : GNU GPL License
// Description : LiftWave Wavelet Library Component
//============================================================================
/*
* Copyright (c) 2012 Rafat Hussain
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*/
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

	} else if (name == "db5") {

		    ltype="dpdpdp";
			stages=6;
			Kconst=1.2314418287192634;
			
			//Stage 1,2,3,4,5,6
			
			double d1[]={-0.265145};
			double p1[]={0.247729,-0.878163};
			double d2[]={0.534125,0.241421};
			double p2[]={0.1985336258386243,-0.6332784120192370};
			double d3[]={-0.0877884834474499,0.0137333394082371};
			double p3[]={-0.0315951369981596};
			
			lcoeff.insert(lcoeff.begin(),p3,p3+1);
			lcoeff.insert(lcoeff.begin(),d3,d3+2);
			lcoeff.insert(lcoeff.begin(),p2,p2+2);
			lcoeff.insert(lcoeff.begin(),d2,d2+2);
			lcoeff.insert(lcoeff.begin(),p1,p1+2);
			lcoeff.insert(lcoeff.begin(),d1,d1+1);
			
			int pl[]={1,0,2,0,2,1,2,1,2,-1,1,2};
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

	} else if (name == "bior2.6") {

		    ltype="dp";
			stages=2;
			Kconst=0.707107;
			
			//Stage 1,2
			
			double d1[]={0.5,0.5};
			double p1[]={-.00976563,0.0761719,-0.316406,-0.316406,0.0761719,-.00976563};
			
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+6);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,6,2};
			plen.assign (pl,pl + sizeof(pl)/sizeof(int));

	} else if (name == "bior2.8") {

		    ltype="dp";
			stages=2;
			Kconst=0.707107;
			
			//Stage 1,2
			
			double d1[]={0.5,0.5};
			double p1[]={0.00213623,-0.0204468,0.0953979,-0.327087,-0.327087,0.0953979,
							-0.0204468,0.00213623};
			
			
			
			lcoeff.insert(lcoeff.begin(),p1,p1+8);
			lcoeff.insert(lcoeff.begin(),d1,d1+2);
			
			int pl[]={2,1,8,3};
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
			vector<double> poly_coeff;
			for (int j=0; j < plen[2*i];j++) {
				cout << lcoeff[total+j] << " ";
				poly_coeff.push_back(lcoeff[total+j]);
				t2++;
			}
			total=total+t2;
			cout << endl;
			cout << "Laurent Polynomial : ";
			Laurent<double> polydisp;
			polydisp.setPoly(poly_coeff,plen[2*i+1]);
			polydisp.dispPoly();
			cout << endl;
	}
	cout << "--------------------------" << endl;
	cout << " K : " << Kconst <<endl;
}	

void addLift(string &c,vector<double> addcoeff, int mp) {
	ltype=ltype+c;
	stages=ltype.size();
	
	int len_add=addcoeff.size();
	plen.push_back(len_add);
	plen.push_back(mp);
	
	lcoeff.insert(lcoeff.end(),addcoeff.begin(),addcoeff.end());
}

virtual ~liftscheme() {
	
}	
	
};

#endif // LIFT_H
