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

struct indxl{

list<int> elements;
int Size=0;

};

class Matrix{
public:

	unordered_map<int, indxl> row;
	unordered_map<int, indxl> column;
	unordered_map<pair<int, int>, double, pair_hash> M;

	Matrix(){

	}

	Matrix(unordered_map<pair<int, int>, double, pair_hash> M, unordered_map<int, indxl> row, unordered_map<int, indxl> column){
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
	for (auto x : column[i].elements)
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
			row[x].elements.push_back(j);
			row[x].elements.remove(i);
		}

	}

	c1.second = i;
	c2.second = j;
	for (auto x : column[j].elements)
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
			row[x].elements.push_back(i);
			row[x].elements.remove(j);
		}

	}

	s.clear();

	swap(column[i], column[j]);


}


void Matrix::rowSwap(int i, int j){

	pair<int, int> c1, c2;
	unordered_set<int> s;

	c1.first = i;
	c2.first = j;
	for (auto x : row[i].elements)
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
			column[x].elements.push_back(j);
			column[x].elements.remove(i);
		}

	}

	c1.first = i;
	c2.first = j;
	for (auto x : row[j].elements)
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
			column[x].elements.push_back(i);
			column[x].elements.remove(j);
		}

	}

	s.clear();

	swap(row[i], row[j]);


}

void Matrix::print(){

	cout << row.size() << ' ' << column.size()<< ' ' << M.size() << endl;
	for (auto x : M)
		cout << x.first.first << ' ' << x.first.second << ' ' << x.second << endl;


}

void Matrix::fprint(string filename){
	ofstream out(filename);
	out << row.size()<< ' ' << column.size() << ' ' << M.size() << endl;
	for (auto x : M)
		out << x.first.first << ' ' << x.first.second << ' ' << x.second << endl;
	out.close();
}


class LUSOL{
private:
	Matrix U;
	Matrix L;
	int n;
	int m;

	vector<double> b;
	vector<double> ans;
	deque<pair<int, int>> P;
	deque<pair<int, int>> Q;

	void Factorize(Matrix A);
	void pivot(int i);
	void solve(Matrix A,vector<double> a);

public:

	LUSOL(){


	}

	LUSOL(Matrix A, vector<double> a, int N, int M){

        n=N;
        m=M;
		Factorize(A);
		solve(A, a);


	}

	~LUSOL(){

	}


};


void LUSOL::pivot(int i){

	pair<int, int> ix;
	ix.first=i;
	int M2=(U.column[i].Size - 1)*(U.row[i].Size - 1), m2=i;

	for (auto x : U.row[i].elements)
	{

		ix.second= x;
		if (x >=i){

			if (M2 > (U.column[x].Size - 1)*(U.row[i].Size - 1))
			{
				M2 = (U.column[x].Size - 1)*(U.row[i].Size - 1);
				m2 = x;
			} else if ((U.column[x].Size - 1)*(U.row[i].Size - 1)==M2) {
                int d=0;
                for(auto y: U.column[x].elements)
                    if(y-i>d){
                        M2=(U.column[x].Size - 1)*(U.row[i].Size - 1);
                        m2 = x;
                        d=y-i;
                    }
			}
		}


	}


	if(i!=m2) {

        U.colSwap(i,m2);
        Q.push_back({i,m2});
	}


	pair<int, int> xi({ i, i });
	int m1 = i, M1 = (U.column[i].Size - 1)*(U.row[i].Size - 1);

	xi.second = i;

	for (auto x : U.column[i].elements)
	{
		xi.first = x;
		if (x >=i){

			if (M1 > (U.column[i].Size - 1)*(U.row[x].Size - 1))
			{
				M1= (U.column[i].Size - 1)*(U.row[x].Size - 1);
				m1 = x;
			} else if ((U.column[i].Size - 1)*(U.row[x].Size - 1)==M1) {
                int d=0;
                for(auto y: U.row[x].elements)
                    if(y-i>d){
                        M1= (U.column[i].Size - 1)*(U.row[x].Size - 1);
                        m1 = x;
                        d=y-i;
                    }
			}

		}

	}
	if(i!=m1){

        U.rowSwap(i, m1);
        swap(b[i], b[m1]);
        L.colSwap(i, m1);
        L.rowSwap(i, m1);
        P.push_back({i,m1});
	}


}



void LUSOL::Factorize(Matrix A){

    U=A;
    pair<int, int> ij;

    for (int i = 1; i <= n; i++)
    {
        ij.first = ij.second = i;
        L.M[ij] = 1;
        L.row[ij.first].elements.push_back(ij.second);
        L.row[ij.first].Size++;
        L.column[ij.second].elements.push_back(ij.first);
        L.column[ij.second].Size++;
    }

	pair<int, int> ji, ii;

	auto t1 = high_resolution_clock::now();
	for (int i = 1; i <= n-1; i++)
	{
		ii.first = ii.second = i;
		ji.second = i;

		list<int> column = U.column[i].elements;

		pivot(i);

		for (auto j : column)
		{
			ji.first = j;
			if (j>i &&  U.M.find(ji) != U.M.end())
			{

				double p = U.M[ji];

				L.M[ji] = U.M[ji] / U.M[ii];
				L.column[i].elements.push_back(j);
				L.column[i].Size++;
				L.row[j].elements.push_back(i);
				L.row[j].Size++;

				list<int> row = U.row[i].elements;
				for (auto k : row)
				if (k >= i){

					pair<int, int> jk({ j, k });
					pair<int, int> ik({ i, k });

					if (U.M.find(jk) == U.M.end()){
						U.M[jk] = -p*U.M[ik] / U.M[ii];
						U.row[j].elements.push_back(k);
						U.row[j].Size++;
						U.column[k].elements.push_back(j);
						U.column[k].Size++;
					}
					else if (fabs(U.M[jk] - p*U.M[ik] / U.M[ii]) < 0.1) {
						U.M.erase(jk);
						U.row[j].elements.remove(k);
						U.row[j].Size--;
						U.column[k].elements.remove(j);
						U.column[k].Size--;
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

	L.fprint("L.txt");
    U.fprint("U.txt");

}


void LUSOL::solve(Matrix A, vector<double> a){

    b.resize(n+1);
    ans.resize(n+1);

    for (int i = 1; i <= n; i++)
    {
        for (auto x : A.row[i].elements)
        {
            pair<int, int> ix({ i, x });
            b[i] += A.M[ix]*a[x];
        }
    }



	pair<int, int> ij;

	vector<double> y(n+1);

	auto t1 = high_resolution_clock::now();
	for (int i = 1; i <=n; i++)
	{

		y[i] = b[i];
		ij.first = i;
		for (auto x : L.row[i].elements)
		if (x<i){
			ij.second = x;
			y[i] -= L.M[ij] * y[x];
		}


	}



	for (int i = n; i >= 1; i--)
	{

		ans[i] = y[i];
		ij.first = i;
		for (auto x : U.row[i].elements)
		if (x>i){
			ij.second = x;
			ans[i] -= U.M[ij] * ans[x];
		}
		ij.second = i;
		ans[i] /= U.M[ij];
	}

	while(!Q.empty()){

        ij=Q.back();
        Q.pop_back();
        swap(ans[ij.first], ans[ij.second]);
	}

	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> ms_double = t2 - t1;
	ofstream out("SolveTime.txt");
	out << ms_double.count() << "ms";
	out.close();


	double error = 0;
	for (int i = 1; i <= n; i++)
		error += (ans[i] - a[i])*(ans[i] - a[i]);
	error = sqrtf(error);

	ofstream out1("solution.txt");
	out1 << "Error: " << error<<"\n\n";
	for (int i = 1; i <= n; i++)
		out1 <<a[i]<<'\t' << ans[i] << '\n';
	out1.close();

}


int main()
{


	int nonzrel, n, m;
	double value;
	pair<int, int> ij;
	Matrix A;


	ifstream in("1138_bus.mtx");
	if(!in.is_open()) {cout<<"The file doesn't exist"; return 0;}

	in >> n >> m >> nonzrel;
    vector<double> ans(n+1);

	for (int i = 1; i <= nonzrel; i++)
	{
		in >> ij.first >> ij.second >> value;
		A.M[ij] = value;
		A.row[ij.first].elements.push_back(ij.second);
		A.row[ij.first].Size++;
		A.column[ij.second].elements.push_back(ij.first);
		A.column[ij.second].Size++;
	}


    for(int i=1; i<=n; i++)
        ans[i]=rand()%100;


	LUSOL lusol(A, ans, n, m);


	return 0;
}
