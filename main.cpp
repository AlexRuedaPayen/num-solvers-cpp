#include <iostream>
#include <fstream>
#include "Matrix.hpp"
#include "Vector.hpp"
#include <random>

using namespace std;

CSR Laplacien_CSR(int n) {
	CSR G(n,n);
	for (int i=2;i<n;i++) {G(i,i-1)=-1;G(i,i+1)=-1;G(i,i)=2;}
	G(1,1)=2;G(1,2)=-1;G(n,n)=2;G(n,n-1)=-1;
	return G;}

 SparseMatrix Laplacien_SparseMatrix(int n) {
		SparseMatrix G(n,n);
		for (int i=2;i<n;i++) {G(i,i-1)=-1;G(i,i+1)=-1;G(i,i)=2;}
		G(1,1)=2;G(1,2)=-1;G(n,n)=2;G(n,n-1)=-1;
		return G;}


int main () {
mt19937 gen(time(NULL));
normal_distribution<double> N(0,1);
	
fstream Files("MatrixFiles");
string Matrix;
vector<double> vec;

JSON Output1("./output/CGvsPCG.json");
JSON Output2("./output/LUvsMR.json");

std::cout<<"\nData processing"<<"\n";


//Output1.start("LaplaceCSR");
//Output2.start("LaplaceCSR");

//std::cout<<"LaplaceCSR"<<"\n";

CSR B(Laplacien_CSR(100));

vec.clear();

for (int i=0;i<B.dim().first;i++) {
vec.push_back(N(gen));
}

//B.Question1(vec,Output1);
//B.Question2(vec,Output2);

//Output1.end();
//Output2.end();

Output1.start("LaplaceSparse");
Output2.start("LaplaceSparse");

std::cout<<"LaplaceSparse"<<"\n";

SparseMatrix C(Laplacien_SparseMatrix(100));

vec.clear();

for (int i=0;i<C.dim().first;i++) {
vec.push_back(N(gen));
}
C.Question1(vec,Output1);
C.Question2(vec,Output2);

Output1.end();
Output2.end();


while(Files>>Matrix) {
fstream Data(Matrix);
DenseMatrix A(Data);
std::cout<<Matrix<<"\n";

vec.clear();

for (int i=0;i<A.dim().first;i++) {
vec.push_back(N(gen));
}

Output1.start(Matrix);
Output2.start(Matrix);

A.Question2(vec,Output2);
A.Question1(vec,Output1);

Output1.end();
Output2.end();
}


Output1.close();
Output2.close();

Files.close();



return 0;}

