#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "Laurent.h"
#include "filter.h"
#include "lift.h"

using namespace std;

int main()
{
	
	string name="db3";
	liftscheme lft(name);
	vector<double> coeff;
	vector<int> lenv;
	string lat;
	double K;
	lft.getScheme(coeff,lenv,lat,K);
	cout << lat << " " << K << " " << coeff.size() << " " << lenv.size() << endl;
	

	return 0;
}