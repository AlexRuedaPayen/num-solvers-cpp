#include "Matrix.hpp"

SparseMatrix::SparseMatrix(int a,int b) :nr(a),nc(b),Lower(nullptr),Upper(nullptr) {}
SparseMatrix::SparseMatrix(const SparseMatrix & A) : nr(A.nr),nc(A.nc),val(A.val),Lower(A.Lower),Upper(A.Upper) {}
SparseMatrix & SparseMatrix::operator=(const SparseMatrix & A) {
		nr=A.nr;nc=A.nc;val=A.val;Lower=A.Lower;Upper=A.Upper;
	return (*this);
}
double SparseMatrix::operator() (int i,int j) const {
	auto A=val.find(std::make_pair(i,j));
	if(A==val.end()) {return 0;}
	else {return A->second;}}

double & SparseMatrix::operator() (int i,int j) {
	return(val[std::make_pair(i,j)]);
}

std::pair<int,int> SparseMatrix::dim() const {return std::make_pair(nr,nc);}

std::ostream & SparseMatrix::operator>>(std::ostream & o) const {
	o<<"\n";
	auto it=val.begin();
	for (int i=1;i<=nr;i++) {
		for (int j=1;j<=nc;j++){
			if ((it->first.first==i) & (it->first.second==j)) {o<<it->second<<"\t";++it;} else{o<<0<<"\t";} 
		} 
		o<<"\n";
	} o<<"\n";
	return o;
}

void SparseMatrix::MvProd(const Vector & A,Vector & B) const {
	for (auto it=val.begin();it!=val.end();++it) {
		B[(it->first).first-1]+=A[(it->first).second-1]*(it->second);
	}
} 

void SparseMatrix::solve_triangle_inf(const Vector & B, Vector & A) const {
	auto it=val.begin();
	double s;
	while(it!=val.end()) {
		s=0;
		int a=(it->first).first;
		while (it->first.second!=a) {
			s+=A[(it->first).second-1]*(it->second);
			++it;
		}
		A[a-1]=(B[(it->first).second-1]-s)/(it->second);
		++it;
	}
}

void SparseMatrix::solve_triangle_sup(const Vector & B, Vector & A) const {
	auto it=val.end();
	double s;
	do{	--it;
		s=0;
		int a=it->first.first;
		while (it->first.second!=a) {
			s+=A[(it->first).second-1]*(it->second);
			--it;
		}
		A[a-1]=(B[(it->first).second-1]-s)/(it->second);
	}
	while (it!=val.begin());}

void SparseMatrix::LU1() {
	Lower=new SparseMatrix(*this);
	Upper=new SparseMatrix(nr,nc);
	auto it1=Upper->val.begin();
	for (int i=1;i<=nr;i++) {
		it1=(Upper->val.insert(it1,std::make_pair(std::make_pair(i,i),1)));
		auto it5=it1;
		auto it2=Lower->val.find(std::make_pair(i,i));
		auto it3=it2;
		++it1;++it2;
		while (it2->first.first==i) {
			it1=(Upper->val.insert(it1,std::make_pair(std::make_pair(i,it2->first.second),it2->second/it3->second)));
			it2=Lower->val.erase(it2);
			++it1;
		}
		++it5;
		for (int j=i+1;j<=nr;j++) {
			while (it5->first.first==i) {
				auto it4=Lower->val.find(std::make_pair(j,i));
				if(it4!=val.end()) {(*Lower)(j,it5->first.second)-=it4->second*it5->second;++it5;}
				else {break;}
			}
		}
	}
}

void SparseMatrix::LUSolve(const Vector & A, Vector & B) {
	if(Lower==nullptr) {LU1();}
	Vector C(B);
	Lower->solve_triangle_inf(A,C);
	Upper->solve_triangle_sup(C,B);
}

SparseMatrix::SparseMatrix(std::istream & is) {
	is>>nr>>nc;
	auto it=val.begin();
	double a;
	for (int i=1;i<=nr;i++) {
		for (int j=1;j<=nc;j++) {
			is>>a;
			if (a!=0) {
				val.insert(it,std::make_pair(std::make_pair(i,j),a));
				++it;
			}
		}
	}
}


int invV1(const std::vector<int> & v, const int i) {
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
}
int invV2(const std::vector<int> & v, const int i) {
	int n=v.size();
	int min=0;
	int max=n-1;
	int a=(n-1)/2;
	if (v[n-1]<i) {return n;}
	if (v[0]>=i) {return 0;}
	else {
		while (i>v[a] | v[a-1]>=i) {
			if (i>=v[a]) {a+=(max-a+1)/2;min=a;} else {a-=(a-min)/2;max=a;}
		}
	return a;
	}
}/*
SparseMatrix SparseMatrix::operator()(const std::vector<int> & Ir,const std::vector<int> & Ic) const {
	SparseMatrix A(Ir.size(),Ic.size());
	int i=0;
	int j=0;
	for (auto it=val.begin();it!=val.end();++it) {
		if (it->first.first>Ir[i]) {
			i=invV2(Ir,it->first.first);
		}
		if (i==Ir.size()) {break;}
		if (it->first.first==Ir[i]) {
			if (it->first.second>Ic[j]) {j=invV1(Ic,it->first.second);}
			if (it->first.second==Ic[j]) {A(i+1,j+1)=it->second;}
		}
	}
	return A; 
}
*/
SparseMatrix SparseMatrix::operator()(const std::vector<int> & Ir,const std::vector<int> & Ic) const {
	SparseMatrix A(Ir.size(),Ic.size());
	for (int i=0;i<Ir.size();i++) {
		for (int j=0;j<Ic.size();j++) {
			double a=this->operator()(Ir[i],Ic[j]);
			if (a!=0) {A(i+1,j+1)=a;}
		}
	}
	return(A);
}

SparseMatrix::~SparseMatrix() {
	delete Lower;
	delete Upper;
}


SparseMatrix operator*(const SparseMatrix & A,const SparseMatrix & B) {
	SparseMatrix C(A.dim().first,B.dim().second);
	for (auto it1=A.val.begin();it1!=A.val.end();++it1) {
		for (auto it2=B.val.begin();it2!=B.val.end();++it2) {
			if (it1->first.second==it2->first.first) {
				C(it1->first.first,it2->first.second)+=(it1->second)*(it2->second);
			}
		}
	}
	return C;
 }



