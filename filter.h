#ifndef FILTER_H_
#define FILTER_H_
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "Laurent.h"

using namespace std;

void biorfilt(string ,Laurent<double> &,Laurent<double> &,Laurent<double> &,
                            Laurent<double> &) ;

void orthfilt(string ,Laurent<double> &,Laurent<double> &,Laurent<double> &,
                            Laurent<double> &);
void lpoly(string ,Laurent<double> &,Laurent<double> &,Laurent<double> &,
                            Laurent<double> &);							


#endif /* FILTER_H_ */
