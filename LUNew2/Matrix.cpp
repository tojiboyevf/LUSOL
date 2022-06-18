#include "Matrix.h"

Matrix::Matrix()
{

}

Matrix::Matrix(std::vector<std::map<int, double>> M, std::vector<std::set<int>> column){
    this->column = column;
    this->M = M;
}


Matrix::~Matrix()
{

}

void Matrix::colSwap(int i, int j){

	std::set<int> s;

	for (auto x : column[i])
	{
		s.insert(x);

		if (M[x].find(j) != M[x].end()) {
			std::swap(M[x][i], M[x][j]);
		}
		else if (M[x].find(j) == M[x].end()){
			M[x][j] = M[x][i];
			M[x].erase(i);
		}

	}

	for (auto x : column[j])
	if (s.count(x) == 0)
	{

		if (M[x].find(i) != M[x].end()) {
			std::swap(M[x][j], M[x][i]);
		}
		else if (M[x].find(i) == M[x].end()){
			M[x][i] = M[x][j];
			M[x].erase(j);
		}

	}

	s.clear();

	std::swap(column[i], column[j]);


}


void Matrix::rowSwap(int i, int j){


	std::set<int> s;
    std::map<int, double> row=M[i];
	for (auto x : row)
	{
		int y= x.first;
		s.insert(y);

		if (M[j].find(y) != M[j].end()) {
			std::swap(M[j][y], M[i][y]);
		}
		else if (M[j].find(y)== M[j].end()){
			M[j][y] = M[i][y];
			M[i].erase(y);
			column[y].insert(j);
			column[y].erase(i);
		}

	}

	row=M[j];
	for (auto x : row)
	if (s.count(x.first) == 0)
	{
        int y=x.first;

		if (M[i].find(y) != M[i].end()) {
			std::swap(M[i][y], M[j][y]);
		}
		else if (M[i].find(y) == M[i].end()){
			M[i][y] = M[j][y];
			M[j].erase(y);
			column[y].insert(i);
			column[y].erase(j);
		}

	}
    row.clear();
	s.clear();


}

void Matrix::print(){

    int s=0;
	for (int i=1; i<M.size(); i++)
        s+=M[i].size();

	std::cout << M.size()-1 << ' ' << column.size()-1 << ' ' << M.size() << std::endl;
	for (int i=1; i<M.size(); i++)
    {
        for(auto x: M[i])
            std::cout <<i<< ' ' << x.first << ' ' << x.second << std::endl;
    }



}

void Matrix::fprint(std::string filename){

	int s=0;
	for (int i=1; i<M.size(); i++)
        s+=M[i].size();

    std::ofstream out(filename);
    out << M.size()-1 << ' ' << column.size()-1 << ' ' << s << std::endl;
	for (int i=1; i<M.size(); i++)
    {
        for(auto x: M[i])
            out <<i<< ' ' << x.first << ' ' << x.second << std::endl;
    }
	out.close();
}
