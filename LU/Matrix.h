#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <fstream>


class Matrix
{
public:

    std::vector<std::set<int>> column;
	std::vector<std::map<int, double>> M;

	Matrix();

	Matrix(std::vector<std::map<int, double>> M, std::vector<std::set<int>> column);

	~Matrix();

	void colSwap(int i, int j);
	void rowSwap(int i, int j);
	void fprint(std::string file);
	void print();
};

#endif // MATRIX_H
