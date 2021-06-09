#include <iostream>
#include"Matrix.h"
#include"LUd.h"

int main()
{
    int nonzrel, n, m,i,j;
	double value;

	Matrix A;


	std::ifstream in("1138.mtx");
	if(!in.is_open()) {std::cout<<"The file doesn't exist"; return 0;}

	in >> n >> m >> nonzrel;

    A.column.resize(m+1);
    A.M.resize(n+1);
	for (int k = 1; k <= nonzrel; k++)
	{
		in >> i>> j>> value;
		A.M[i][j] = value;
		A.column[j].insert(i);
	}


	LUd lu(A, n, m);

	return 0;
}
