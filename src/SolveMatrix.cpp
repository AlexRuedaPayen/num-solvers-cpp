
#include "Matrix.hpp"
#include<numeric>
#include<cmath>
#include<iostream>
#include<ctime>
#include "JSON.hpp"
#include "Block_code.hpp"

std::list<VectorInt> gen_list(int l,int r,int n) {
	std::list<VectorInt> L;
	VectorInt V;
	for (int k=1;k<=l+r;k++){
		V.push_back(k);
	}
	L.push_back(V);
	for (int j=1;j<=((n/l)-2);j++) {
		V.clear();
		for (int k=l*j-r;k<=l*(j+1)+r;k++)
			V.push_back(k);
		L.push_back(V)
		;
	}
	V.clear();
	for (int k=(l*(n/l)-1);k<=n;k++){
		V.push_back(k);
	}
	L.push_back(V);
	return(L);
}



void SolveMatrix::cg(const Vector & A, Vector & B, double tol, int nmax,JSON & json1) const {
	json1.start("cg");
	double Res=norm2(A);
	double epsilon=tol*Res;
	Vector P(A);
	Vector R(A);
	double beta;
	double alpha;
	double Res2;
	Vector Zero(A.size());
	Vector Prod(A.size());
	Vector ResVecIt;
	int n=0;
	ResVecIt.push_back(Res);
	while ((Res>=epsilon) && (n<nmax)) {
		n+=1;
		Prod=Zero;
		MvProd(P,Prod);
		alpha=((Res*Res)/(P*Prod));
		Operator1(B,alpha,P);
		Res2=Res;
		Operator2(R,alpha,Prod);
		Res=sqrt(R*R);
		beta=((Res*Res)/(Res2*Res2));
		Operator3(P,beta,R);
		ResVecIt.push_back(Res);
	}
	json1.add("resvec",ResVecIt);
	json1.end();
}

void SolveMatrix::MinRes(const Vector & A,Vector & B,double tol,int nmax) const {
	Vector R(A);
	double Res=sqrt(R*R);
	double epsilon=tol*Res;
	double alpha;
	int n=0;
	Vector Zero(A.size());
	Vector Prod1(Zero);
	Vector Prod2(Zero);
	while ((Res>=epsilon) && (n<nmax)) {
		Prod1=Zero;
		Prod2=Zero;
		MvProd(R,Prod1);
		alpha=((R*Prod1)/(Prod1*Prod1));
		Operator1(B,alpha,R);
		MvProd(B,Prod2);
		R=A-Prod2;
		Res=sqrt(R*R);
		n+=1;
	}
}


void SolveMatrix::Question1(const Vector & A,JSON & json1)const {
	Vector B(A.size());
	cg(A,B,1e-6,dim().first,json1);
	B=Vector(A.size());
	pcg(A,B,1e-6,dim().first,json1);
	}

void SolveMatrix::Question2(const Vector & A,JSON & json2) {
	Vector B(A.size());
	clock_t t1,t2,t3,t4;
	t1=clock();
	LUSolve(A,B);
	t2=clock();
	Vector C(A.size());
	t4=clock();
	MinRes(A,C,1e-6,1e6);
	t3=clock();
	json2.start("LUSolve");
	json2.add("Time",static_cast<double>(t2-t1) /CLOCKS_PER_SEC);
	json2.end();
	json2.start("MinRes");
	json2.add("Time",static_cast<double>(t3-t4) /CLOCKS_PER_SEC);
	json2.end();
}