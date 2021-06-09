#include "LUd.h"

LUd::LUd()
{

}


LUd::LUd(Matrix A, int n, int m)
{
   this->U = A;
   this->n=n;
   this->m=m;
   L.column.resize(n+1);
   L.M.resize(n+1);
   this->Factorize();
   this->solve();
}

LUd::~LUd()
{

}


void LUd::pivot(int i){

    int m;

    auto e=c.begin();
    std::pair<int,int> x=*e;

    while(U.M[i].find(x.second)==U.M[i].end()){
        e++;
        x=*e;
    }

    c.erase(x);
    if(i!=x.second){
        m=x.second;
        c.erase({U.column[i].size(),i});
        c.insert({U.column[i].size(),m});
        U.colSwap(i,m);
        Q.push_back({i,m});

    }


    auto t=r.begin();
    x=*t;
    while(U.column[x.second].count(i)==0){
        t++;
        x=*t;
    }
    r.erase(x);

    if(i!=x.second){
        m=x.second;
        r.erase({U.M[i].size(),i});
        r.insert({U.M[i].size(),m});
        U.rowSwap(i, m);
        std::swap(b[i], b[m]);
        L.rowSwap(i,m);
        P.push_back({i,m});
    }


}



void LUd::Factorize(){


	a.resize(n+1);
	ans.resize(n+1);
	b.resize(n+1);

	for(int i=1; i<=n; i++)
        a[i]=rand()%100;

	for (int i = 1; i <= n; i++)
	{
		for (auto x : U.M[i])
		{
			b[i] += x.second*a[x.first];
		}
	}



    for(int i=1; i<=n; i++)
    {
        c.insert({U.column[i].size(), i});
        r.insert({U.M[i].size(), i});
    }

	auto t1 = high_resolution_clock::now();
	for (int i = 1; i <n; i++)
	{

		std::set<int> column = U.column[i];

		pivot(i);

		for (auto j : column)
		{
			if (j>i &&  U.M[j].find(i) != U.M[j].end())
			{

				double p = U.M[j][i];

				L.M[j][i] = U.M[j][i] / U.M[i][i];


				std::map<int, double> row = U.M[i];
				for (auto K : row)
				if (K.first >= i){

                    int k=K.first;

					if (U.M[j].find(k) == U.M[j].end()){

                        if(r.count({U.M[j].size(),j})){
                            r.erase({U.M[j].size(),j});
                            r.insert({U.M[j].size()+1,j});
                        }
						U.M[j][k] = -p*U.M[i][k] / U.M[i][i];

                        if(c.count({U.column[k].size(),k})){
                            c.erase({U.column[k].size(),k});
                            c.insert({U.column[k].size()+1,k});
                        }
						U.column[k].insert(j);

					}
					else if (fabs(U.M[j][k] - p*U.M[i][k] / U.M[i][i]) < 0.001) {

                        if(r.count({U.M[j].size(),j})){
                            r.erase({U.M[j].size(),j});
                            r.insert({U.M[j].size()-1,j});
                        }

						U.M[j].erase(k);

                        if(c.count({U.column[k].size(),k})){
                            c.erase({U.column[k].size(),k});
                            c.insert({U.column[k].size()-1,k});
                        }

						U.column[k].erase(j);

					}
					else if (U.M[j].find(k) != U.M[j].end()){
						U.M[j][k] -= p*U.M[i][k] / U.M[i][i];
					}

				}
				row.clear();

			}
		}
		column.clear();
	}
	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> ms_double = t2 - t1;
	std::ofstream out("FactorTime.txt");
	out << ms_double.count() << "ms";
	out.close();


}


void LUd::solve(){


    for (int i = 1; i <=n; i++)
		{
			L.M[i][i] = 1;
		}

    L.fprint("L.txt");
	U.fprint("U.txt");

	std::vector<double> y(n+1);


	auto t1 = high_resolution_clock::now();
	for (int i = 1; i <=n; i++)
	{

		y[i] = b[i];
		for (auto x : L.M[i])
		if (x.first<i){
			y[i] -= L.M[i][x.first] * y[x.first];
		}


	}



	for (int i = n; i >= 1; i--)
	{

		ans[i] = y[i];
		for (auto x : U.M[i])
		if (x.first>i){
			ans[i] -= U.M[i][x.first] * ans[x.first];
		}
		ans[i] /= U.M[i][i];
	}

    std::pair<int, int> ij;

	while(!Q.empty()){

        ij=Q.back();
        Q.pop_back();
        std::swap(ans[ij.first], ans[ij.second]);
	}



	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> ms_double = t2 - t1;
	std::ofstream out("SolveTime.txt");
	out << ms_double.count() << "ms";
	out.close();


	double error = 0;
	for (int i = 1; i <= n; i++)
		error += (ans[i] - a[i])*(ans[i] - a[i]);
	error = sqrtf(error);

	std::ofstream out1("solution.txt");
	out1 << "Error: " << error<<"\n\n";

	for (int i = 1; i <= n; i++)
		out1 <<a[i]<<'\t' <<ans[i] << '\n';

	out1.close();


}




