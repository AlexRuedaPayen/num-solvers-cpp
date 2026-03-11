#include<numeric>
#include<cmath>
#include<vector>
#include<map>
#include<list>
#include<ostream>

#ifndef VECTOR_
#define VECTOR_

typedef std::pair<int,int> Index;
struct comp {bool operator()(const Index & a,const Index & b) const {if (a.first==b.first){return(a.second<b.second);} else {return(a.first<b.first);}}};
typedef std::map<Index,double,comp> Identity;
typedef std::pair<int,double> ColVal;
typedef std::list<ColVal>::iterator Iter;
typedef std::vector<int> VectorInt;
typedef std::vector<double> Vector;

void Operator1(Vector &, double,const Vector &);
void Operator2(Vector &, double,const Vector &);
void Operator3(Vector &, double,const Vector &);
double norm2(const Vector &);
double operator*(const Vector &, const Vector &);
Vector operator+(const Vector &, const Vector &);
Vector operator-(const Vector &, const Vector &);
std::ostream & operator<<(std::ostream &, const Vector &);
std::ostream & operator<<(std::ostream &, const VectorInt &);

#endif
