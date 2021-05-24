#include <iostream>
#include <unordered_map>
#include <list>
#include <vector>
#include <deque>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


using namespace std;

struct pair_hash
{
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2> &pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

class Matrix{
public:

	unordered_map<int, list<int>> row;
	unordered_map<int, list<int>> column;
	unordered_map<pair<int, int>, double, pair_hash> M;

	Matrix(){

	}

	Matrix(unordered_map<pair<int, int>, double, pair_hash> M, unordered_map<int, list<int>> row, unordered_map<int, list<int>> column){
		this->row = row;
		this->column = column;
		this->M = M;

	}

	~Matrix(){

	}

	void colSwap(int i, int j);
	void rowSwap(int i, int j);
	void fprint(string file);
	void print();


};

void Matrix::colSwap(int i, int j){

	pair<int, int> c1, c2;
	unordered_set<int> s;

	c1.second = i;
	c2.second = j;
	for (auto x : column[i])
	{
		s.insert(x);
		c1.first = x;
		c2.first = x;
		if (M.find(c2) != M.end()) {
			swap(M[c2], M[c1]);
		}
		else if (M.find(c2) == M.end()){
			M[c2] = M[c1];
			M.erase(c1);
			row[x].push_back(j);
			row[x].remove(i);
		}

	}

	c1.second = i;
	c2.second = j;
	for (auto x : column[j])
	if (s.count(x) == 0)
	{
		c1.first = x;
		c2.first = x;
		if (M.find(c1) != M.end()) {
			swap(M[c2], M[c1]);
		}
		else if (M.find(c1) == M.end()){
			M[c1] = M[c2];
			M.erase(c2);
			row[x].push_back(i);
			row[x].remove(j);
		}

	}

	s.clear();

	list<int> c;
	swap(column[i], column[j]);

	c.clear();

}


void Matrix::rowSwap(int i, int j){

	pair<int, int> c1, c2;
	unordered_set<int> s;

	c1.first = i;
	c2.first = j;
	for (auto x : row[i])
	{
		s.insert(x);
		c1.second = x;
		c2.second = x;
		if (M.find(c2) != M.end()) {
			swap(M[c2], M[c1]);
		}
		else if (M.find(c2) == M.end()){
			M[c2] = M[c1];
			M.erase(c1);
			column[x].push_back(j);
			column[x].remove(i);
		}

	}

	c1.first = i;
	c2.first = j;
	for (auto x : row[j])
	if (s.count(x) == 0)
	{

		c1.second = x;
		c2.second = x;
		if (M.find(c1) != M.end()) {
			swap(M[c2], M[c1]);
		}
		else if (M.find(c1) == M.end()){
			M[c1] = M[c2];
			M.erase(c2);
			column[x].push_back(i);
			column[x].remove(j);
		}

	}

	s.clear();

	list<int> c;
	swap(row[i], row[j]);

	c.clear();

}

void Matrix::print(){

	cout << row.size() << ' ' << column.size() << ' ' << M.size() << endl;
	for (auto x : M)
		cout << x.first.first << ' ' << x.first.second << ' ' << x.second << endl;


}

void Matrix::fprint(string filename){
	ofstream out(filename);
	out << row.size() << ' ' << column.size() << ' ' << M.size() << endl;
	for (auto x : M)
		out << x.first.first << ' ' << x.first.second << ' ' << x.second << endl;
	out.close();
}


class LUSOL{
private:
	void Factorize();
	void pivot(int i);
	void solve();

public:
	Matrix U;
	Matrix L;
	vector<double> X;
	vector<double> ans;
	deque<pair<int, int>> P;

	LUSOL(){


	}

	LUSOL(Matrix A, vector<double> b, int n){
		this->U = A;

		pair<int, int> ij;

		for (int i = 1; i <= n; i++)
		{
			ij.first = ij.second = i;
			this->L.M[ij] = 1;
			this->L.row[ij.first].push_back(ij.second);
			this->L.column[ij.second].push_back(ij.first);
		}

		this->X = b;
		this->ans.resize(n+1);
		this->Factorize();
		this->solve();

	}

	~LUSOL(){

	}


};


void LUSOL::pivot(int i){


	pair<int, int> xi({ i, i });
	int m = 0, M = -1;

	if (U.M.find(xi) == U.M.end())
	{
		bool b = true;
		for (auto x : U.column[i])
		{
			xi.first = x;
			if (x >= i && U.M.find(xi) != U.M.end() && U.M[xi] != 0){
				m = x;
				b = false;
				M = (U.column[i].size() - 1)*(U.row[m].size() - 1);

			}
			if (!b) break;
		}


	}

	if (M == -1){ m = i; M = (U.column[i].size() - 1)*(U.row[i].size() - 1); }
	xi.second = i;

	for (auto x : U.column[i])
	{
		xi.first = x;
		if (x >i && U.M.find(xi) != U.M.end() && U.M[xi] != 0){

			if (M > (U.column[i].size() - 1)*(U.row[x].size() - 1))
			{
				M = (U.column[i].size() - 1)*(U.row[x].size() - 1);
				m = x;
			}
		}

	}

	U.rowSwap(i, m);
	swap(X[i], X[m]);
	L.colSwap(i, m);
	L.rowSwap(i, m);
    P.push_back({i,m});

}



void LUSOL::Factorize(){

	int dim = U.row.size();
	pair<int, int> ji, ii;

	auto t1 = high_resolution_clock::now();
	for (int i = 1; i <= dim - 1; i++)
	{
		ii.first = ii.second = i;
		ji.second = i;

		list<int> column = U.column[i];

		pivot(i);

		for (auto j : column)
		{
			ji.first = j;
			if (j>i &&  U.M.find(ji) != U.M.end())
			{

				double p = U.M[ji];

				L.M[ji] = U.M[ji] / U.M[ii];
				L.column[i].push_back(j);
				L.row[j].push_back(i);

				list<int> row = U.row[i];
				for (auto k : row)
				if (k >= i){

					pair<int, int> jk({ j, k });
					pair<int, int> ik({ i, k });

					if (U.M.find(jk) == U.M.end()){
						U.M[jk] = -p*U.M[ik] / U.M[ii];
						U.row[j].push_back(k);
						U.column[k].push_back(j);
					}
					else if (fabs(U.M[jk] - p*U.M[ik] / U.M[ii]) < 0.1) {
						U.M.erase(jk);
						U.row[j].remove(k);
						U.column[k].remove(j);
					}
					else if (U.M.find(jk) != U.M.end()){
						U.M[jk] -= p*U.M[ik] / U.M[ii];
					}

				}
				row.clear();

			}
		}
		column.clear();
	}
	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> ms_double = t2 - t1;
	ofstream out("FactorTime.txt");
	out << ms_double.count() << "ms";
	out.close();

}


void LUSOL::solve(){

	int dim = X.size();
	pair<int, int> ij;

	vector<double> y(dim);

	auto t1 = high_resolution_clock::now();
	for (int i = 1; i <dim; i++)
	{

		y[i] = X[i];
		ij.first = i;
		for (auto x : L.row[i])
		if (x<i){
			ij.second = x;
			y[i] -= L.M[ij] * y[x];
		}


	}



	for (int i = dim-1; i >= 1; i--)
	{

		ans[i] = y[i];
		ij.first = i;
		for (auto x : U.row[i])
		if (x>i){
			ij.second = x;
			ans[i] -= U.M[ij] * ans[x];
		}
		ij.second = i;
		ans[i] /= U.M[ij];
	}
	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> ms_double = t2 - t1;
	ofstream out("SolveTime.txt");
	out << ms_double.count() << "ms";
	out.close();

}



int main()
{


	int nonzrel, n, m;
	double value;
	pair<int, int> ij;
	Matrix A;


	/*string filenameaddress;
	cin>>filenameaddress;*/

	ifstream in("5000x5000.txt");
	if(!in.is_open()) {cout<<"The file doesn't exist"; return 0;}

	in >> n >> m;

	in >> nonzrel;

	for (int i = 1; i <= nonzrel; i++)
	{
		in >> ij.first >> ij.second >> value;
		A.M[ij] = value;
		A.row[ij.first].push_back(ij.second);
		A.column[ij.second].push_back(ij.first);
	}

    vector<double> ans(n+1);
	vector<double> b(n+1);


	/*for (int i = 1; i <= n; i++)
		cin >> b[i];*/

	/*for(int i=1; i<=n; i++)
        cin>>ans[i];*/

    for(int i=1; i<=n; i++)
        ans[i]=rand()%100;

	for (int i = 1; i <= n; i++)
	{
		for (auto x : A.row[i])
		{
			pair<int, int> ix({ i, x });
			b[i] += A.M[ix]*ans[x];
		}
	}

	LUSOL lusol(A, b, n);
	lusol.L.fprint("L.txt");
	lusol.U.fprint("U.txt");

	double error = 0;
	for (int i = 1; i <= n; i++)
		error += (lusol.ans[i] - ans[i])*(lusol.ans[i] - ans[i]);
	error = sqrtf(error);

	ofstream out("solution.txt");
	out << "Error: " << error<<"\n\n";

	for (int i = 1; i <= n; i++)
		out <<ans[i]<<'\t' << lusol.ans[i] << '\n';

	out.close();
	return 0;
}
