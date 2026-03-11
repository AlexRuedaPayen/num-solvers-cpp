#include "Vector.hpp"

void Operator1(Vector & A, double alpha, const Vector & B) {
for (int i=0;i<A.size();i++) {A[i]+=(alpha*B[i]);}
}

void Operator2(Vector & A, double alpha, const Vector & B) {
for (int i=0;i<A.size();i++) {A[i]-=(alpha*B[i]);}
}

void Operator3(Vector & A, double alpha,const Vector & B) {
for (int i=0;i<A.size();i++) {A[i]=(B[i]+(alpha*A[i]));}
}

double norm2(const Vector & A) {
	return(sqrt(std::accumulate(A.begin(),A.end(),0,[](double s,double a){return(s+a*a);})));
}

double operator*(const Vector & A, const Vector & B) {
double s=0;
for (int i=0;i<A.size();i++) {s+=(A[i]*B[i]);}
return s;
}

Vector operator+(const Vector & A, const Vector & B) {
Vector C(B); 
for (int i=0;i<A.size();i++) {C[i]=(A[i]+B[i]);}
return C;
}
Vector operator-(const Vector & A, const Vector & B) {
Vector C(B); 
for (int i=0;i<A.size();i++) {C[i]=(A[i]-B[i]);}
return C;
}

std::ostream & operator<<(std::ostream & o, const Vector & A) {
	for (int i=0; i<A.size();i++) {
		o<<A[i]<<"\t";
	}
	o<<"\n";
	return o;
}

std::ostream & operator<<(std::ostream & o, const VectorInt & A) {
	for (int i=0; i<A.size();i++) {
		o<<A[i]<<"\t";
	}
	o<<"\n";
	return o;
}