#include "Matrix.hpp"
#include <cmath>
#include <iostream>

DenseMatrix::DenseMatrix(int a, int b) : nr(a),nc(b),val(new double[a*b]),LU(nullptr){
for (int i=0;i<nr*nc;i++)
{*(val+i)=0;}
}

DenseMatrix::DenseMatrix(const DenseMatrix & A) : nr(A.nr),nc(A.nc),val(new double[nr*nc]),LU(A.LU){
for (int i=0;i<nr*nc;i++)
{*(val+i)=*((A.val)+i);}
}

DenseMatrix::DenseMatrix(std::istream & i) : val(nullptr),LU(nullptr){
std::string a,b;
i>>a>>nr>>nc>>b;
val=new double[nr*nc];
for (int j=0;j<nr*nc;j++) {
	i>>*(val+j);
}
}


DenseMatrix & DenseMatrix::operator=(const DenseMatrix & A) {
	if (this!=&A) {
		delete[] val;
		delete[] LU;
		nr=A.nr;
		nc=A.nc;
		val = new double[nr*nc];
		for (int i=0;i<nr*nc;i++) {
			*(val+i)=*((A.val)+i);
		}
		if (A.LU!=nullptr) {
		LU = new DenseMatrix(nr,nc);
			for (int i=0;i<nr*nc;i++) {
				*(LU+i)=*((A.LU)+i);
			}
		}}
		return(*this);
	}

double DenseMatrix::operator() (int a,int b) const {
	return(*(val+(nc*(a-1)+(b-1))));
}


double & DenseMatrix::operator() (int a,int b){
	delete LU;
	return(*(val+(nc*(a-1)+(b-1))));
}

std::ostream & DenseMatrix::operator>>(std::ostream & o) const {
	o<<"\n";
	for (int i=1;i<=nr;i++) {
		for (int j=1;j<=nc;j++) {
			o<<(*(val+(nc*(i-1)+(j-1))))<<"\t";
		}
		o<<"\n";
	}
	o<<"\n";
	return o;
}

void DenseMatrix::MvProd (const Vector & A, Vector & B) const {
	for (int i=1;i<=nc;i++) {
		for (int j=1;j<=nr;j++) {
			B[i-1]+=(*(val+(nc*(i-1)+(j-1))))*A[j-1];
		}
	}
}


std::pair<int,int> DenseMatrix::dim() const {
return(std::make_pair(nr,nc))
;}

void DenseMatrix::solve_triangle_inf(const Vector & Y,Vector & X) const {
for (int i=0; i<nc; i++) {double s=0;
for (int j=0; j<i; j++) {s+=(*(val+((nc*i)+j)))*X[j];}
X[i]=(Y[i]-s)/(*(val+((nc*i)+i)));}}

void DenseMatrix::solve_triangle_sup(const Vector & Y,Vector & X) const {
for (int i=nc-1; i>=0; i--) {double s=0;
for (int j=nc-1; j>i; j--) {s+=(*(val+((nc*i)+j)))*X[j];}
X[i]=(Y[i]-s)/(*(val+((nc*i)+i)));
}}

void DenseMatrix::LU1(void) {
LU=new DenseMatrix(*this);
for (int i=1;i<=nr;i++) {for (int j=1+i;j<=nc;j++) 
{(*LU)(i,j)=(*LU)(i,j)/(*LU)(i,i);}
{for (int j=i+1;j<=nr;j++) {
	for (int k=i+1;k<=nc;k++) {(*LU)(j,k)=(*LU)(j,k)-((*LU)(j,i)*(*LU)(i,k));}}}}}

void DenseMatrix::LUSolve(const Vector & B,Vector & A)
{if (nc==nr) {
	if (LU==nullptr) {LU1();}
	Vector C(B),buffer;
	LU->solve_triangle_inf(B,C);
	for (int i=1;i<=nr;i++) {
	buffer.push_back((*LU)(i,i));
	(*LU)(i,i)=1;}
	 LU->solve_triangle_sup(C,A);
	for (int i=1;i<=nr;i++) {
	(*LU)(i,i)=buffer[i-1];};
	}}

DenseMatrix DenseMatrix::operator()(const std::vector<int> & Ir,const std::vector<int> & Ic) const {
	DenseMatrix A(Ir.size(),Ic.size());
	for (int i=0;i<Ir.size();i++) {
		for (int j=0;j<Ic.size();j++) {
			A(i+1,j+1)=*(val+((nc*(Ir[i]-1))+(Ic[j]-1)));
		}
	}
	return(A);
}


DenseMatrix::~DenseMatrix() {
	delete[] val;
	delete LU;
}

DenseMatrix operator*(const DenseMatrix & A,const DenseMatrix & B) {
	DenseMatrix C(A.dim().first,B.dim().second);
	for (int i=1;i<=A.dim().second;i++) {
		for (int j=1;j<=B.dim().second;j++) {
			for (int k=1;j<=A.dim().first;k++) {
				C(k,j)+=A(k,i)*B(i,j);
			}
		}
	}
	return C;
} 


