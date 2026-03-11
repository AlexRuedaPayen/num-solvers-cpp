#include<list>
#include<vector>
#include<utility>
#include<istream>
#include<map>
#include<vector>
#include<random>
#include "JSON.hpp"
#include "Vector.hpp"

#ifndef MATRIX_
#define MATRIX_

std::list<VectorInt> gen_list(int l,int r,int n) ;

template<class Matrix>
class Block;

template<class Matrix>
class BlockMatrix;

class SolveMatrix {
	private:
		virtual void pcg(const Vector &,Vector &,double,int,JSON &) const=0;
		void cg(const Vector &,Vector &,double,int,JSON &) const;
	public :
		void MinRes(const Vector &,Vector &,double,int) const;
		void Question1(const Vector &,JSON &)const;
		void Question2(const Vector &,JSON &);
		virtual void MvProd(const Vector &,Vector &) const=0; 
		virtual std::pair<int,int> dim() const=0;
		virtual void LUSolve(const Vector &,Vector &)=0;
		//virtual ~SolveMatrix()=0;
};

class DenseMatrix : public SolveMatrix {
	private:
		int nr,nc;
		double * val;
		void LU1();
		DenseMatrix * LU;
		void pcg(const Vector &,Vector &,double,int,JSON &) const;
		double & operator() (int,int);
	public:
		friend Block<DenseMatrix>;
		friend BlockMatrix<DenseMatrix>;
		DenseMatrix(int,int);
		void solve_triangle_inf(const Vector &,Vector &) const;
 	    void solve_triangle_sup(const Vector &,Vector &) const;
 	    void LUSolve(const Vector &,Vector &);
		DenseMatrix(const DenseMatrix &);
		DenseMatrix(std::istream &);
		DenseMatrix & operator=(const DenseMatrix &);
		~DenseMatrix();
		double operator() (int,int) const;
		std::pair<int,int> dim() const;
		std::ostream & operator>>(std::ostream & o) const;
		void MvProd (const Vector &,Vector &) const;
		DenseMatrix operator()(const std::vector<int> &,const std::vector<int> &) const;
		friend DenseMatrix operator*(const DenseMatrix &,const DenseMatrix &); 
};

class SparseMatrix : public SolveMatrix {
 	private :
 		int nr,nc;
 		void pcg(const Vector &,Vector &,double,int,JSON &) const;
 		Identity val;
 		SparseMatrix * Lower;
 		SparseMatrix * Upper;
 		void LU1();
 		double & operator() (int,int);
 	public :
 		friend Block<SparseMatrix>;
 		friend BlockMatrix<SparseMatrix>;
		SparseMatrix(int,int);
		void solve_triangle_inf(const Vector &, Vector &) const;
 		void solve_triangle_sup(const Vector &, Vector &) const;
		SparseMatrix(const SparseMatrix &);
		SparseMatrix(std::istream &);
		~SparseMatrix();
		SparseMatrix & operator=(const SparseMatrix &);
		double operator() (int,int) const;
		std::ostream & operator>>(std::ostream & o) const;
		void MvProd (const Vector &,Vector &) const;
		std::pair<int,int> dim() const;
		void LUSolve(const Vector &,Vector &);
		SparseMatrix operator()(const std::vector<int> &,const std::vector<int> &) const;
		friend SparseMatrix operator*(const SparseMatrix &,const SparseMatrix &); 
		friend SparseMatrix Laplacien_SparseMatrix(int n);
};

class CSR : public SolveMatrix {
	private :
		int nr,nc;
		std::list<ColVal> colval;
		CSR * Upper;
 		CSR * Lower;
 		void LU1();
		std::vector<Iter> row;
		void pcg(const Vector &,Vector &,double,int,JSON &) const;
		double & operator() (int,int);
	public:
		friend Block<CSR>;
		friend BlockMatrix<CSR>;
		CSR(int,int);
		void solve_triangle_inf(const Vector &, Vector &) const;
 		void solve_triangle_sup(const Vector &, Vector &) const;
 		~CSR();
		CSR(const CSR &);
		CSR & operator=(const CSR &);
		double operator()(int,int) const;
		std::ostream & operator>>(std::ostream & o) const;
		void MvProd(const Vector &,Vector &) const;
		std::pair<int,int> dim() const;
		void LUSolve(const Vector &,Vector &);
		CSR operator()(const std::vector<int> &,const std::vector<int> &) const;
		friend CSR operator*(const CSR &,const CSR &); 
		friend CSR Laplacien_CSR(int n);
};

int invV(const std::vector<int> &, const int);
int invV1(const std::vector<int> &,const int);
int invV2(const std::vector<int> &,const int);
#endif
