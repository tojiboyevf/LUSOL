# LUSOL
This project factorizes given matrix into Lower(L) and Upper(U) triangular matrices. To check how it is working correctly it generates random dense vector and solves the system of linear equations
# Structure

## tests
``tests`` folder contains different test cases. Matrices are given in the coordinate text format. In coordinate text file format the first line lists three integers: the number of rows `m`, columns `n`, and nonzeros `nz` in the matrix. The nonzero matrix elements are then listed, one per line, by specifying row index `i`, column index `j`, and the value `a(i,j)`, in that order. For example,

```
     m       m       nz
     i1      j1      val1
     i2      j2      val2
     i3      j3      val3
     .       .        .
     .       .        .
     .       .        .
     inz     jnz     valnz
    
```
File ``5000x5000.txt`` also contains a matrix in the coordinate text format.
## main.cpp
There are two main classes: ``Matrix`` and ``LUSOL`` in ``main.cpp``. <br>
### Matrix
``Matrix`` has three main attributes and two main methods: <br>
``unordered_map<int, list<int>> row`` - row[i] contains indexes of nonzero elements in i th row <br>
``unordered_map<int, list<int>> column`` - column[j] contains indexes of nonzero elements in j th column <br>
``unordered_map<pair<int, int>, double, pair_hash> M`` - for key value it is better to use ``pair``  and mapped value is in ``double`` type <br>
Method ``void colSwap(int i, int j)`` swaps i th column and j th column<br>
Method ``void rowSwap(int i, int j)`` swaps i th row and j th row<br>

### LUSOL
``LUSOL`` has five main attributes and three main methods: <br>
``Matrix U`` is Upper triangular matrix <br>
``Matrix L`` is Lower triangular matrix <br>
``vector<double> b`` is column vector <br>
``vector<double> ans`` is solution of the system of linear equation <br>
``deque<pair<int, int>> P`` contains permutation matrices  <br>
Method ``void Factorize()`` factorizes given matrix ``A`` into matrices ``L`` and ``U`` <br>
Method ``void pivot(int i)`` during factorization chooses such element ``A[ij]`` with minimum ``(A.column[i].size() - 1)*(A.row[m].size() - 1)`` for preserving sparsity of matrix<br>
Method ``void solve()`` solves the system of linear equations<br>
