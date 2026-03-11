#include<vector>
#include<list>
#include "Matrix.hpp"

#ifndef BLOCK_
#define BLOCK_

template <class Matrix>
class Block {
	private:
		int nr,nc;
		VectorInt Ir,Ic;
		Matrix mat;
		void MvProdInv(const Vector &, Vector &);
	public:
		friend class BlockMatrix<Matrix>;
		Block(int,int,const VectorInt &,const VectorInt &);
		Block(const Block<Matrix> &);
		Block & operator=(const Block<Matrix> &);
		void MvProd (const Vector &, Vector &) const;
	};

template <class Matrix>
class BlockMatrix {
	private:
		int nr,nc;
		std::list<Block<Matrix>> val;
		void MvProdInv(const Vector &, Vector &);
	public:
		friend class Block<Matrix>;
		friend class DenseMatrix;
		friend class CSR;
		friend class SparseMatrix;
		BlockMatrix(int,int);
		std::pair<int,int> dim() const {return(std::make_pair(nr,nc));}
		BlockMatrix(const BlockMatrix<Matrix> &);
		BlockMatrix & operator=(const BlockMatrix<Matrix> &);
		BlockMatrix & operator+=(const Block<Matrix> &);
		void MvProd(const Vector &,Vector &) const;
		void Extract(const Matrix &,const std::list<VectorInt> &);
};
#endif