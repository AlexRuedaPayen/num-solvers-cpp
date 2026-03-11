#include "Matrix.hpp"
#include<iostream>

CSR::CSR(int a,int b) : nr(a),nc(b),colval({std::make_pair(nc+1,0)}),Upper(nullptr),Lower(nullptr) {
	row=std::vector<Iter>(nr+1,colval.begin());}

CSR::CSR(const CSR & A) : nr(A.nr),nc(A.nc),colval(A.colval),Upper(A.Upper),Lower(A.Lower),row(A.row) {}

CSR & CSR::operator=(const CSR & A) {
	row=A.row;colval=A.colval;Lower=A.Lower;Upper=A.Upper;
		return *this;}

double CSR::operator() (int i,int j) const {
	Iter a=row[i-1];
	while (a!=row[i]) {
	if (a->first==j) {return a->second;} 
	a++;}
	return 0;
}

CSR::~CSR() {
	delete Lower;
	delete Upper;
}

double & CSR::operator() (int i,int j) {
	Iter a=row[i-1];
	if (a->first>j) {
		row[i-1]=colval.insert(a,std::make_pair(j,0));
		return(row[i-1]->second);
	}
	else { 
		while (a!=row[i]) {
		if (a->first==j) {
			return a->second;
		} 
		if (a->first>j) {
			return (colval.insert(a,std::make_pair(j,0))->second);
		}
			++a;
	}
	return(colval.insert(a,std::make_pair(j,0))->second);
	}
}

std::pair<int,int> CSR::dim() const {
	return std::make_pair(nr,nc);}

std::ostream & CSR::operator>>(std::ostream & o) const {
	o<<"\n";
	Iter it;
	for (int i=1;i<=nr;i++) {
		it=row[i-1];
		for (int j=1;j<=nc;j++) {
			if (((it->first)==j) & (it!=row[i])) {o<<it->second<<"\t";++it;}
			else {o<<0<<"\t";}
		}
		o<<"\n";
	}
return o;}

void CSR::MvProd(const Vector & A,Vector & B) const {
	for (int i=1;i<=nr;i++) {
		for (Iter it=row[i-1];it!=row[i];++it) {
			B[i-1]+=(it->second)*A[it->first -1];
		}
}}

void CSR::solve_triangle_inf(const Vector & A,Vector & B) const {
	for (int i=1;i<=nr;i++) {
		double s=0;
		Iter it=row[i-1];
		while(it->first!=i) {
			s+=B[it->first -1]*(it->second);
			++it;
		}
		B[i-1]=(A[it->first -1]-s)/(it->second);
	}
}

void CSR::solve_triangle_sup(const Vector & A,Vector & B) const {
		for (int i=nr;i>=1;i--) {
		double s=0;
		Iter it=row[i];
		--it;
		while(it->first!=i) {
			s+=B[it->first -1]*(it->second);
		--it;
		}
		B[i-1]=(A[it->first-1]-s)/(it->second);
	}
}


void CSR::LU1() {
	Lower=new CSR(*this);
	Upper=new CSR(nr,nc);
	for (int i=1;i<=nr;i++) {
		auto it1=Lower->row[i-1];
		while (it1->first!=i) {++it1;}
		auto it2=Upper->row[i-1];
		while (it2->first!=i) {++it2;}
		auto it3=it1;
		it2=Upper->colval.insert(it2,std::make_pair(i,1));
		auto it4=it2;
		Upper->row[i-1]=it2;
		++it1;++it2;
		while (it1->first>i) {
			it2=Upper->colval.insert(it2,std::make_pair(it1->first,it1->second/it3->second));
			it1=Lower->colval.erase(it1);
			++it2;
		}
		++it4;
		for (int j=i+1;j<=nr;j++) {
			auto it5=Lower->row[j];
			while (it5!=Lower->row[j+1]) {++it5;if(it5->first>=i) {break;}}
			if (it5->first==i) {
			while (it4->first>i) {
			(*Upper)(j,it4->first)-=it5->second*it4->second;++it4;
			}}
		}
	}
}

void CSR::LUSolve(const Vector & A, Vector & B) {
	if(Lower==nullptr & Upper==nullptr) {LU1();}
	Vector C(B);
	Lower->solve_triangle_inf(A,C);
	Upper->solve_triangle_sup(C,B);
}

int invV(const std::vector<int> & v, const int i) {
	int n=v.size();
	int min=0;
	int max=n-1;
	int a=(n-1)/2;
	if (v[n-1]<i) {return 0;}
	if (v[0]>=i) {return 0;}
	else {
		while (i>v[a] | v[a-1]>=i) {
			if (i>=v[a]) {a+=(max-a+1)/2;min=a;} else {a-=(a-min)/2;max=a;}
		}
	return a;
	}
}/*

CSR CSR::operator()(const std::vector<int> & Ir,const std::vector<int> & Ic) const {
	CSR A(Ir.size(),Ic.size());
	for (int i=0;i<Ir.size();i++) {
		for (auto it=row[Ir[i]-1];it!=row[Ir[i]];++it) {
			int j=0;
			if (it->first>Ic[j]) {j=invV(Ic,it->first); if (j==0) {break;}}
			if (it->first==Ic[j]) {A(i+1,j+1)=it->second;}
		}
	}
	return A;
}
*/

CSR CSR::operator()(const std::vector<int> & Ir,const std::vector<int> & Ic) const {
	CSR A(Ir.size(),Ic.size());
	for (int i=0;i<Ir.size();i++) {
		for (int j=0;j<Ic.size();j++) {
			double a=this->operator()(Ir[i],Ic[j]);
			if (a!=0) {A(i+1,j+1)=a;}
		}
	}
	return(A);
}
