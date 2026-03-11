#include "Block.hpp"
#include<numeric>
#include<algorithm>

#ifndef BLOCK_CODE_
#define BLOCK_CODE_

template<class Matrix>
Block<Matrix>::Block(int a,int b,const VectorInt & I,const VectorInt & J) :  nr(a), nc(b), Ir(I), Ic(J), mat(I.size(),J.size()) {}


template<class Matrix>
Block<Matrix>::Block(const Block<Matrix> & A) : nr(A.nr),nc(A.nc),Ir(A.Ir),Ic(A.Ic),mat(A.mat) {}

template<class Matrix>
Block<Matrix> & Block<Matrix>::operator=(const Block<Matrix> & A) {
	nr=A.nr;nc=A.nc;Ir=A.Ir;Ic=A.Ic;mat=A.mat;return(*this);
}

template<>
void Block<DenseMatrix>::MvProd(const Vector & A,Vector & B) const {
	for (int i=0;i<Ir.size();i++) {
		for (int j=0;j<Ic.size();j++) {
			B[Ir[i]-1]+=mat(i+1,j+1)*A[Ic[j]-1];
		}
	}
}

template<>
void Block<SparseMatrix>::MvProd(const Vector & A,Vector & B) const {
	for (auto it=mat.val.begin();it!=mat.val.end();++it) {
		B[Ir[it->first.first-1]-1]+=(it->second)*A[Ic[it->first.second-1]-1];
	}
}

template<>
void Block<CSR>::MvProd(const Vector & A,Vector & B) const {
	for (int i=1;i<=Ir.size();i++) {
		for (auto it=mat.row[i-1];it!=mat.row[i];++it) {
			B[Ir[i-1]-1]+=(it->second)*A[Ic[it->first-1]-1];
		}
	}
}

template<class Matrix>
BlockMatrix<Matrix>::BlockMatrix(int i,int j) : nr(i),nc(j) {}

template<class Matrix>
BlockMatrix<Matrix>::BlockMatrix(const BlockMatrix & A) : nr(A.nr),nc(A.nc),val(A.val) {}

template<class Matrix>
BlockMatrix<Matrix> & BlockMatrix<Matrix>::operator=(const BlockMatrix<Matrix> & A) {
	nr=A.nr;nc=A.nc;val=A.val;return(*this);
}


template<class Matrix>
BlockMatrix<Matrix> & BlockMatrix<Matrix>::operator+=(const Block<Matrix> & A) {
	val.push_back(A);return(*this);
}

template<class Matrix>
void BlockMatrix<Matrix>::MvProd(const Vector & A,Vector & B) const {
	B=std::accumulate(val.begin(),val.end(),B,[& A](Vector & B,Block<Matrix> beta){
		Vector i(A.size());
		beta.MvProd(A,i);
		return B+i;}); }

template<class Matrix>
void BlockMatrix<Matrix>::Extract(const Matrix & A,const std::list<VectorInt> & List) {
	nr=A.nr;
	nc=A.nc;
	for (auto it=List.begin();it!=List.end();++it) {
		Matrix B=A(*it,*it);
		Block<Matrix> Beta(A.nr,A.nc,*it,*it);
		Beta.mat=B;
		this->operator+=(Beta);
	}
}

template<>
void Block<DenseMatrix>::MvProdInv(const Vector & A,Vector & B) {
	Vector C(Ic.size());
	Vector D(Ic.size());
	for (int i=0;i<Ic.size();i++) {C[i]=A[Ic[i]-1];}
	mat.LUSolve(C,D);
	for (int i=0;i<Ir.size();i++) {B[Ir[i]-1]=D[i];}
}


template<>
void Block<SparseMatrix>::MvProdInv(const Vector & A,Vector & B) {
	Vector C(Ic.size());
	Vector D(Ic.size());
	for (int i=0;i<Ic.size();i++) {C[i]=A[Ic[i]-1];}
	mat.MinRes(C,D,1e-6,1e4);
	for (int i=0;i<Ir.size();i++) {B[Ir[i]-1]=D[i];}
}


template<>
void Block<CSR>::MvProdInv(const Vector & A,Vector & B) {
	Vector C(Ic.size());
	Vector D(Ic.size());
	for (int i=0;i<Ic.size();i++) {C[i]=A[Ic[i]-1];}
	mat.MinRes(C,D,1e-6,1e4);
	for (int i=0;i<Ir.size();i++) {B[Ir[i]-1]=D[i];}
}


template<class Matrix>
void BlockMatrix<Matrix>::MvProdInv(const Vector & A,Vector & B) {
	B=std::accumulate(val.begin(),val.end(),B,[& A](Vector & B,Block<Matrix> beta){
		Vector i(A.size());
		beta.MvProdInv(A,i);
		return B+i;});
}

void CSR::pcg(const Vector & A,Vector & B,double tol,int nmax,JSON & json1) const {
	json1.start("pcg");
	std::mt19937 gen(time(NULL));
	std::uniform_int_distribution<int> C36PO(2,floor(dim().first/2));
	//int l=C36PO(gen);
	//std::uniform_int_distribution<int> R2D2(1,l-1);
	//int r=R2D2(gen);
	int l=log(dim().first)/log(2),r=l-1;
	BlockMatrix<CSR> Bl(dim().first,dim().first);
	std::list<VectorInt> L=gen_list(l,r,dim().first);
	Bl.Extract(*this,gen_list(l,r,dim().first));
	std::vector<int> Lis({l,r,dim().first});
	double Res=sqrt(A*A);
	double epsilon=tol*Res;
	Vector Zero(B.size());
	Vector PA(Zero);
	Vector P(Zero);
	Vector R(A);
	Bl.MvProdInv(R,P);
	Vector Z(P);
	double beta;
	double gamma;
	double alpha;
	Vector ResVecIt;
	int n=0;
	ResVecIt.push_back(Res);
	while ((Res>=epsilon) && (n<nmax)) {
		n+=1;
		PA=Zero;
		MvProd(P,PA);
		gamma=R*Z;
		alpha=gamma/(P*PA);
		Operator1(B,alpha,P);
		Operator2(R,alpha,PA);
		Z=Zero;
		Bl.MvProdInv(R,Z);
		Res=sqrt(R*R);
		beta=((R*Z)/gamma);
		Operator3(P,beta,Z);
		ResVecIt.push_back(Res);
	}
	json1.add("parameters",Lis);
	json1.add("resvec",ResVecIt);
	json1.end();
}


void SparseMatrix::pcg(const Vector & A,Vector & B,double tol,int nmax,JSON & json1) const {
	json1.start("pcg");
	int l=log(dim().first)/log(2),r=l-1;
	BlockMatrix<SparseMatrix> Bl(dim().first,dim().first);
	std::list<VectorInt> L=gen_list(l,r,dim().first);
	Bl.Extract(*this,gen_list(l,r,dim().first));
	std::vector<int> Lis({l,r,dim().first});
	double Res=sqrt(A*A);
	double epsilon=tol*Res;
	Vector Zero(B.size());
	Vector PA(Zero);
	Vector P(Zero);
	Vector R(A);
	Bl.MvProdInv(R,P);
	Vector Z(P);
	double beta;
	double gamma;
	double alpha;
	Vector ResVecIt;
	int n=0;
	ResVecIt.push_back(Res);
	while ((Res>=epsilon) && (n<nmax)) {
		n+=1;
		PA=Zero;
		MvProd(P,PA);
		gamma=R*Z;
		alpha=gamma/(P*PA);
		Operator1(B,alpha,P);
		Operator2(R,alpha,PA);
		Z=Zero;
		Bl.MvProdInv(R,Z);
		Res=sqrt(R*R);
		beta=((R*Z)/gamma);
		Operator3(P,beta,Z);
		ResVecIt.push_back(Res);
	}
	json1.add("parameters",Lis);
	json1.add("resvec",ResVecIt);
	json1.end();
}

void DenseMatrix::pcg(const Vector & A,Vector & B,double tol,int nmax,JSON & json1) const {
	json1.start("pcg");
	int l=log(dim().first)/log(2);int r=l-1;
	BlockMatrix<DenseMatrix> Bl(dim().first,dim().first);
	std::list<VectorInt> L=gen_list(l,r,dim().first);
	Bl.Extract(*this,gen_list(l,r,dim().first));
	std::vector<int> Lis({l,r,dim().first});
	double Res=sqrt(A*A);
	double epsilon=tol*Res;
	Vector Zero(B.size());
	Vector PA(Zero);
	Vector P(Zero);
	Vector R(A);
	Bl.MvProdInv(R,P);
	Vector Z(P);
	double beta;
	double gamma;
	double alpha;
	Vector ResVecIt;
	int n=0;
	ResVecIt.push_back(Res);
	while ((Res>=epsilon) && (n<nmax)) {
		n+=1;
		PA=Zero;
		MvProd(P,PA);
		gamma=R*Z;
		alpha=gamma/(P*PA);
		Operator1(B,alpha,P);
		Operator2(R,alpha,PA);
		Z=Zero;
		Bl.MvProdInv(R,Z);
		Res=sqrt(R*R);
		beta=((R*Z)/gamma);
		Operator3(P,beta,Z);
		ResVecIt.push_back(Res);
	}
	json1.add("parameters",Lis);
	json1.add("resvec",ResVecIt);
	json1.end();
}
#endif
