#if !defined(__MATRIXNXM__H__)
#define __MATRIXNXM__H__

#include "svd.h"
#include <geometry.H>
#include <string>
#include <functional>
#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include "GA/gl_util.h"

class MatrixNxM
{
private:
	void FreeMemory(double **&d, int r, int c)
	{
		if(d != nullptr)
		{
			for(int i = 0 ; i < rows ; ++i )
				delete [] d[i];
			delete [] d;
			d = nullptr;
		}
	}

	void AllocateMemory(double **&d, int r, int c)
	{
        d = new double*[r];
		for(int i = 0 ; i < r ; ++i )
			d[i] = new double[c];
	}

public:
    double **data;
    int rows;
    int cols;

	const double*	operator [] ( int i ) const { return data[ i ]; }
	double*			operator [] ( int i ) { return data[ i ]; }

	operator double* () { return (double*)data; }
	operator const double* () const { return (const double*)data; }

    MatrixNxM()
	{
        this->rows = 0;
        this->cols = 0;
		this->data = nullptr;
	}

    MatrixNxM(int rows, int cols)
    {
        this->rows = rows;
        this->cols = cols;
        AllocateMemory(data, rows, cols);
    }

    // Copy constructor.
	MatrixNxM(const MatrixNxM& matrixNxM)
    {
        this->rows = matrixNxM.rows;
        this->cols = matrixNxM.cols;
        AllocateMemory(data, rows, cols);
		
        for( int i = 0 ; i < rows ; ++i )
            for( int j = 0 ; j < cols ; ++j )
				this->data[i][j] = matrixNxM.data[i][j];		
    }

	// Move constructor.
	MatrixNxM( MatrixNxM&& matrixNxM) //move cons
	{
        this->rows = matrixNxM.rows;
        this->cols = matrixNxM.cols;
		this->data = matrixNxM.data;
		matrixNxM.data = nullptr;
	}
	
	// Copy Assignment operator.
	MatrixNxM operator = (const MatrixNxM& matrixNxM)
	{		
		if(this->rows != matrixNxM.rows || this->cols != matrixNxM.cols)
		{
			FreeMemory(data, rows, cols);
			this->rows = matrixNxM.rows;
			this->cols = matrixNxM.cols;
			AllocateMemory(data, rows, cols);
		}

        for( int i = 0 ; i < rows ; ++i )
            for( int j = 0 ; j < cols ; ++j )
				this->data[i][j] = matrixNxM.data[i][j];

		return *this;
	}

	// Move assignment operator.
	MatrixNxM& operator=(MatrixNxM&& matrixNxM)
	{
		if (this != &matrixNxM)
		{
			FreeMemory(data, rows, cols);
			this->rows = matrixNxM.rows;
			this->cols = matrixNxM.cols;
			this->data = matrixNxM.data;
			matrixNxM.data = nullptr;
		}
		return *this;
	}

    ~MatrixNxM()
	{
		FreeMemory(data, rows, cols);
	}

    void ZeroMatrix()
    {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
				data[i][j] = 0.0;
    }

    void Identity()
    {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
            {
                if (i == j)
                    data[i][j] = 1.0;
                else
                    data[i][j] = 0.0;
            }
    }

    MatrixNxM CopyMatrix()
    {
        MatrixNxM c(rows, cols);

        for( int i = 0 ; i < rows ; ++i )
            for( int j = 0 ; j < cols ; ++j )
				c.data[i][j] = data[i][j];

        return c;
    }

	MatrixNxM Transpose()
	{
		MatrixNxM c(cols, rows);

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				c.data[j][i] = data[i][j];
		return c;
	}

    double Determinant()
    {
		if( cols == 3 && rows == 3 )
		{
			return( data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
					data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
					data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]) );
		}
		else if( cols == 2 && rows == 2 )
		{
			return data[0][0] * data[1][1] - data[0][1] * data[1][0];
		}

        int N = this->rows;

        double *a = new double[N * N];
        int *pivot = new int [N];

		auto IND = [N](int i, int j) { return i + j * N; };

        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                a[IND(i,j)] = data[i][j];
            }
        }

        int info = dge_fa(N, a, pivot);

        if (info != 0)
            throw std::string("The factorization failed");

        double det = dge_det(N, a, pivot);

		delete [] a;
		delete [] pivot;

		return det;
    }

	void SVD(MatrixNxM &U, MatrixNxM &W, MatrixNxM &V)
	{
		double **u, **v;
		double *w = new double[cols + 1];

		AllocateMemory(u, rows+1, cols+1);
		AllocateMemory(v, cols+1, cols+1);

		for (int i = 1; i <= this->rows; ++i)
			for (int j = 1; j <= this->cols; ++j)
				u[i][j] = this->data[i-1][j-1];

		svdcmp(u, rows, cols, w, v);

		if(U.rows != rows || U.cols != cols)
			U = MatrixNxM(rows, cols);
		if(V.rows != cols || V.cols != cols)
			V = MatrixNxM(cols, cols);
		if(W.rows != cols || W.cols != cols)
			W = MatrixNxM(cols, cols);

		for (int i = 1; i <= this->rows; ++i)
			for (int j = 1; j <= this->cols; ++j)
				U[i-1][j-1] = u[i][j];

		for (int i = 1; i <= this->cols; ++i)
			for (int j = 1; j <= this->cols; ++j)
				V[i - 1][j - 1] = v[i][j];

		W.ZeroMatrix();

		for (int i = 1; i <= this->cols; ++i)
			W[i - 1][i - 1] = w[i];

		FreeMemory(u, rows+1, cols+1);
		FreeMemory(v, cols+1, cols+1);
		delete [] w;
	}

	MatrixNxM SolveSVD(const MatrixNxM& rhs)
	{
		double **u, **v;
		double *w = new double[cols + 1];
		double *b = new double[rows + 1];
		double *x = new double[cols + 1];

		AllocateMemory(u, rows+1, cols+1);
		AllocateMemory(v, cols+1, cols+1);

		const double illConditionedThreshold = 1e-6;

		for (int i = 1; i <= this->rows; ++i)
			for (int j = 1; j <= this->cols; ++j)
				u[i][j] = this->data[i - 1][j - 1];

		svdcmp(u, rows, cols, w, v);

		for (int i = 1; i <= cols; ++i)
			if (w[i] < illConditionedThreshold)
				w[i] = 0.0;

		MatrixNxM res(rhs.rows, rhs.cols);

		for (int j = 0; j < rhs.cols; ++j)
		{
			for (int i = 0; i < rhs.rows; ++i)
			{
				b[i + 1] = rhs[i][j];
			}

			svbksb(u, w, v, rows, cols, b, x);

			for (int i = 0; i < res.rows; ++i)
			{
				res[i][j] = x[i + 1];
			}
		}

		FreeMemory(u, rows+1, cols+1);
		FreeMemory(v, cols+1, cols+1);
		delete [] w;
		delete [] b;
		delete [] x;

		return res;
	}

    MatrixNxM Solve(const MatrixNxM& rhs)
    {
        int N = this->rows;
        int RHS_NUM = rhs.cols;

        double *a = new double[N * (N + RHS_NUM)];

		auto IND = [N](int i, int j) { return i + j * N; };

        for (int i = 0; i < this->rows; ++i)
        {
            for (int j = 0; j < this->cols; ++j)
            {
                a[IND(i,j)] = this->data[i][j];
            }
            for (int j = 0; j < rhs.cols; ++j)
            {
                a[IND(i, this->cols + j)] = rhs[i][j];
            }
        }

        int solution = r8mat_solve(N, RHS_NUM, a);

        if (solution != 0)
            throw std::string("factorization failed. The solutions could not be computed.");

        MatrixNxM res(rhs.rows, rhs.cols);

        for (int i = 0; i < res.rows; ++i)
        {
            for (int j = 0; j < rhs.cols; ++j)
            {
                res[i][j] = a[IND(i, this->cols + j)];
            }
        }
		
		delete [] a;

        return res;
    }

    MatrixNxM SolveLU(const MatrixNxM& rhs)
    {
        int N = this->rows;

        double *a = new double[N * N];
        double *b = new double[N];
        int	*pivot = new int[N];

		auto IND = [N](int i, int j) { return i + j * N; };

        for (int i = 0; i < this->rows; ++i)
        {
            for (int j = 0; j < this->cols; ++j)
            {
				a[IND(i, j)] = this->data[i][j];
            }
        }

        int info = dge_fa(N, a, pivot);

        if (info != 0)
            throw std::string("The factorization failed");

        MatrixNxM res(rhs.rows, rhs.cols);

        for (int j = 0; j < rhs.cols; ++j)
        {
            for (int i = 0; i < rhs.rows; ++i)
            {
                b[i] = rhs[i][j];
            }

            dge_sl(N, a, pivot, b, 0);

            for (int i = 0; i < res.rows; ++i)
            {
                res[i][j] = b[i];
            }
        }

		delete [] a;
		delete [] b;
		delete [] pivot;

        return res;
    }

	c3ga::rotor ToRotor()
	{
		double trace = data[0][0] + data[1][1] + data[2][2] + 1.0f;
		double qw; // scalar coordinate
		double qx; // coordinate for -e2^e3
		double qy; // coordinate for -e3^e1
		double qz; // coordinate for -e1^e2
		if (trace > 0.00001) {
			double s = 0.5 / (double)sqrt(trace);
			qw = 0.25 / s;
			qw = sqrt(trace) * (0.5);
			qx = (data[2][1] - data[1][2]) * s;
			qy = (data[0][2] - data[2][0]) * s;
			qz = (data[1][0] - data[0][1]) * s;
		}
		else {
			if (data[0][0] > data[1][1] && data[0][0] > data[2][2]) {
				double s = 2.0 * (double)sqrt( 1.0 + data[0][0] - data[1][1] - data[2][2]);
				qx = 0.25 * s;
				qy = (data[0][1] + data[1][0]) / s;
				qz = (data[0][2] + data[2][0]) / s;
				qw = (data[1][2] - data[2][1]) / s;
			}
			else if (data[1][1] > data[2][2]) {
				double s = 2.0 * (double)sqrt( 1.0 + data[1][1] - data[0][0] - data[2][2]);
				qx = (data[0][1] + data[1][0]) / s;
				qy = 0.25 * s;
				qz = (data[1][2] + data[2][1]) / s;
				qw = (data[0][2] - data[2][0]) / s;
			}
			else {
				double s = 2.0 * (double)sqrt( 1.0 + data[2][2] - data[0][0] - data[1][1] );
				qx = (data[0][2] + data[2][0]) / s;
				qy = (data[1][2] + data[2][1]) / s;
				qz = 0.25 * s;
				qw = (data[0][1] - data[1][0]) / s;
			}
		}

		double s = (double) sqrt(qw *qw + qx * qx + qy * qy + qz * qz);

		return c3ga::rotor(c3ga::rotor_scalar_e1e2_e2e3_e3e1, qw / s, -qz / s, -qx / s, -qy / s);
	}

};

inline MatrixNxM operator *(const MatrixNxM& a, const MatrixNxM& b) 
{
    MatrixNxM c(a.rows, b.cols);

	c.ZeroMatrix();

    for (int i = 0; i < a.rows; ++i)
        for (int j = 0; j < b.cols; ++j)
			for( int k = 0 ; k < a.cols ; ++k )
				c[i][j] += a[i][k] * b[k][j];

	return c;
}

inline MatrixNxM operator *(const MatrixNxM& a, double b) 
{
    MatrixNxM c(a.rows, a.cols);

    for ( int i = 0 ; i < a.rows ; ++i)
		for ( int j = 0 ; j < a.cols ; ++j )
			c[i][j] = a[i][j] * b;

	return c;
}

inline MatrixNxM operator *(double b, const MatrixNxM& a) 
{
    return a * b;
}

inline MatrixNxM operator +(const MatrixNxM& a, const MatrixNxM& b) 
{
    MatrixNxM c(a.rows, a.cols);

    for (int i = 0; i < a.rows; ++i)
        for (int j = 0; j < a.cols; ++j)
            c[i][j] = a[i][j] + b[i][j];

    return c;
}


#endif
