#ifndef LUD_H
#define LUD_H

#include <map>
#include <vector>
#include <deque>
#include <set>
#include <string>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <chrono>
#include "Matrix.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;




class LUd
{
private:
	void Factorize();
	void pivot(int i);
	void solve();

public:
	Matrix U;
	Matrix L;
	std::vector<double> b;
	std::vector<double> ans;
	std::vector<double> a;
	std::deque<std::pair<int, int>> P;
	std::deque<std::pair<int, int>> Q;
	int n;
	int m;


	LUd();

	LUd(Matrix A, int n, int m);

	~LUd();
};

#endif // LUD_H
