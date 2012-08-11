#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "Laurent.h"
#include "filter.h"

using namespace std;

void rbiorfilt() {
	
}


int main()
{
	cout << "Lifting Demo" << endl;
	string name="db1";
	Laurent<double> lpd,hpd,lpr,hpr;
	lpoly(name,lpd,hpd,lpr,hpr);
	lpr.dispPoly();
	hpr.dispPoly();
	lpd.dispPoly();
	hpd.dispPoly();
	return 0;
}
