package GiciTransform;

import java.util.Arrays;
import GiciAnalysis.ImageCovariance;
import GiciMatrix.MatrixAlgebra;

public class KarhunenLoeveTransform {
	
	/**
	 * The eigen-decomposition of the covariance matrix.
	 * @param covariance
	 * @return a vector with two matrices (the transform matrix and the diagonal matrix).
	 */
	public static float[][][] generateTransform (float[][] covariance_matrix) {		
		try {
			KarhunenLoeveTransform tk = new KarhunenLoeveTransform();
			return tk.doKLT(toCompleteMirrorFill(covariance_matrix));
		} catch (Exception e) {
			e.printStackTrace();
			throw new Error(e.getMessage());
		}
	}

	public static float[][] toCompleteMirrorFill(final float[][] a) {
		/* Assume it's right aligned */
		int i, width = 0;

		for (i = 0; i < a.length; i++) {
			width = Math.max(width, a[i].length);
		}

		float[][] r = new float[a.length][width];

		for (i = 0; i < a.length; i++) {
			for (int j = 0; j < a[i].length; j++) {
				r[j + width - a[i].length][i] =
					r[i][j + width - a[i].length] = a[i][j];
			}
		}

		return r;
	}
	
	/**
	 * Generates a KLT transform, optimal for the given imageSamples.
	 *
	 * @param imageSamples data to train the transform with.
	 * @return transformation matrix.
	 */
	public static float[][] run(final float[][][] imageSamples) {
		float[][] cv = ImageCovariance.generateCovarianceMatrix(imageSamples, 0);
		return generateTransform(cv)[0];
	}

	/**
	 * Generates a KLT transform, optimal for the given imageSamples.
	 *
	 * @param imageSamples data to train the transform with.
	 * @return transformation matrix.
	 */
	public static float[][][] runWithDiagonal(final float[][][] imageSamples) {
		return generateTransform(ImageCovariance.generateCovarianceMatrix(imageSamples, 0));
	}

	float[][][] doKLT(float[][] covariance_matrix) throws Exception {

		assert (covariance_matrix.length > 0 && covariance_matrix[0].length == covariance_matrix.length);

		final int n = covariance_matrix.length;

		assert(noNaNs(covariance_matrix));
		
		//MatrixAlgebra.printMatrix(covariance_matrix);
		
		gsl_matrix m = new gsl_matrix(covariance_matrix);
		covariance_matrix = null;
		
		gsl_vector eval = new gsl_vector(n);
		gsl_matrix evec = new gsl_matrix(n, n);

		gsl_eigen_symmv_workspace w = new gsl_eigen_symmv_workspace(n);

		gsl_eigen_symmv(m, eval, evec, w);

		gsl_eigen_symmv_sort (eval, evec);

		float[][] matrix = MatrixAlgebra.identityC(n);
		float[][] values = MatrixAlgebra.identityC(n);

		for (int i = 0; i < evec.size1; i++) {
			for (int j = 0; j < evec.size2; j++) {
				matrix[i][j] = (float) gsl_matrix_get(evec, i, j);
			}
		}

		for (int i = 0; i < eval.size; i++) {
			values[i][i] = (float) gsl_vector_get(eval, i);
		}

		float[][][] r = {matrix, values};

		assert(noNaNs(matrix));
		assert(noNaNs(values));
		
		return r;
	}

	boolean noNaNs (float[][] a) {
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[i].length; j++) {
				if (Float.isInfinite(a[i][j]) || Float.isNaN(a[i][j]))
					return false;
			}
		}
		
		return true;
	}
	
	
	/* Code below here is adapted from the GSL - GNU Scientific Library (1.9, 1.13) */

	/* Files from which code has been borrowed:
	 * vector/gsl_vector_double.h
	 * vector/swap_source.c
	 * vector/view_source.c
	 * vector/subvector_source.c
	 * matrix/gsl_matrix_double.h
	 * matrix/swap_source.c
	 * matrix/view_source.c
	 * matrix/init_source.c
	 * matrix/rowcol_source.c
	 * cblas/source_nrm2_r.h
	 * cblas/source_scal_r.h
	 * cblas/source_axpy_r.h
	 * cblas/source_syr2.h
	 * cblas/source_symv.h
	 * cblas/source_dot_r.h
	 * cblas/dnrm2.c
	 * cblas/dscal.c
	 * cblas/daxpy.c
	 * cblas/dsyr2.c
	 * cblas/dsymv.c
	 * cblas/ddot.c
	 * cblas/cblas.h
	 * blas/blas.c
	 * linalg/householder.c
	 * linalg/symmtd.c
	 * eigen/sort.c
	 * eigen/symmv.c
	 * eigen/gsl_eigen.h
	 */

	/*
	 * Author list: Gerard Jungman, Brian Gough, Others?
	 */

	/* Copyright (C) 1996-2009 Gerard Jungman, Brian Gough, 
	 * 
	 * This program is free software; you can redistribute it and/or modify
	 * it under the terms of the GNU General Public License as published by
	 * the Free Software Foundation; either version 3 of the License, or (at
	 * your option) any later version.
	 * 
	 * This program is distributed in the hope that it will be useful, but
	 * WITHOUT ANY WARRANTY; without even the implied warranty of
	 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	 * General Public License for more details.
	 * 
	 * You should have received a copy of the GNU General Public License
	 * along with this program; if not, write to the Free Software
	 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
	 */

	class gsl_eigen_symmv_workspace {
		int size;
		double [] d;
		double [] sd;
		double [] gc;
		double [] gs;

		public gsl_eigen_symmv_workspace(final int n) throws Exception {

			if (n == 0) {
				throw new Exception("matrix dimension must be positive integer");
			}

			d = new double[n];
			sd = new double[n];
			gc = new double[n];
			gs = new double[n];

			size = n;
		}
	};

	class gsl_matrix {
		int size1;
		int size2;
		int tda;
		int start;
		double[] data;

		gsl_matrix() {

		}

		gsl_matrix(final int n, final int m) {
			size1 = n;
			size2 = m;
			tda = size1;
			start = 0;

			data = new double[size1 * size2];
		}

		gsl_matrix(float[][] A) {
			assert (A.length > 0);

			size1 = A.length;
			size2 = A[0].length;
			tda = size1;
			start = 0;

			data = new double[size1 * size2];

			for (int i = 0; i < A.length; i++) {
				assert(A[i].length == A[0].length);

				for (int j = 0; j < A[i].length; j++) {
					data[i * size1 + j] = A[i][j];
				}
			}
		}
	};

	class gsl_vector {
		int size;
		int stride;
		int start;
		double[] data;

		gsl_vector() {
		}

		gsl_vector(final double[] l) {
			stride = 1;
			size = l.length;
			start = 0;

			data = Arrays.copyOf(l, l.length);
		}

		gsl_vector(final int n) {
			data = new double[n];
			size = n;
			stride = 1;
		}
	};

	class gsl_vector_view
	{
		gsl_vector vector;

		gsl_vector_view() {

		}

		gsl_vector_view (double[] base, int n) throws Exception
		{
			if (n == 0)
			{
				throw new Exception("vector length n must be positive integer");
			}

			vector = new gsl_vector();

			vector.data = base;
			vector.size = n;
			vector.stride = 1;
		}
	};

	class gsl_matrix_view {
		gsl_matrix matrix;
	}

	double gsl_matrix_get(final gsl_matrix m, final int i,final int j) throws Exception
	{

		if (i >= m.size1 || i < 0)
		{
			throw new Exception("first int out of range");
		}
		else if (j >= m.size2 || i < 0)
		{
			throw new Exception("second int out of range");
		}

		return m.data[m.start + i * m.tda + j] ;
	} 


	void
	gsl_matrix_set(gsl_matrix m, final int i, final int j, final double x) throws Exception
	{

		if (i >= m.size1 || i < 0)
		{
			throw new Exception("first int out of range");
		}
		else if (j >= m.size2 || i < 0)
		{
			throw new Exception("second int out of range");
		}

		m.data[m.start + i * m.tda + j] = x ;
	}

	void gsl_matrix_set_identity (gsl_matrix m)
	{
		int i, j;
		double [] data = m.data;
		final int p = m.size1;
		final int q = m.size2;
		final int tda = m.tda;

		final double zero = 0;
		final double one = 1;

		for (i = 0; i < p; i++)
		{
			for (j = 0; j < q; j++)
			{
				data[m.start + i * tda + j] = ((i == j) ? one : zero);
			}
		}
	}

	gsl_vector_view gsl_matrix_column (gsl_matrix m, final int j) throws Exception {

		if (j >= m.size2) {
			throw new Exception("column int is out of range");
		}

		gsl_vector v = new gsl_vector();
		gsl_vector_view view = new gsl_vector_view();

		v.data = m.data; // + j
		v.start = m.start + j;
		v.size = m.size1;
		v.stride = m.tda;

		view.vector = v;
		return view;
	}

	void gsl_matrix_swap_columns (gsl_matrix m, final int i, final int j) {
		final int size1 = m.size1;

		if (i != j) {
			int p;

			for (p = 0; p < size1; p++) {
				int n = p * m.tda;

				double tmp = m.data[m.start+n+i] ;
				m.data[m.start+n+i] = m.data[m.start+n+j] ;
				m.data[m.start+n+j] = tmp ;
			}
		}
	}

	gsl_matrix_view gsl_matrix_submatrix (gsl_matrix m, final int i, final int j, final int n1, final int n2) throws Exception {
		if (i >= m.size1) {
			throw new Exception("row index is out of range");
		} else if (j >= m.size2) {
			throw new Exception("column index is out of range");
		} else if (n1 == 0) {
			throw new Exception("first dimension must be non-zero");
		} else if (n2 == 0) {
			throw new Exception("second dimension must be non-zero");
		} else if (i + n1 > m.size1) {
			throw new Exception("first dimension overflows matrix");
		} else if (j + n2 > m.size2) {
			throw new Exception("second dimension overflows matrix");
		}

		gsl_matrix_view view = new gsl_matrix_view(); 
		gsl_matrix s = new gsl_matrix();

		s.data = m.data;
		s.start = m.start + (i * m.tda + j);
		s.size1 = n1;
		s.size2 = n2;
		s.tda = m.tda;

		view.matrix = s;

		return view;
	}

	double
	gsl_vector_get (final gsl_vector v, final int i) throws Exception
	{

		if (i >= v.size || i < 0)
		{
			throw new Exception("int out of range");
		}

		return v.data[v.start + i * v.stride];
	}


	void
	gsl_vector_set (gsl_vector v, final int i, final double x) throws Exception
	{

		if (i >= v.size || i < 0)
		{
			throw new Exception("int out of range");
		}

		v.data[v.start + i * v.stride] = x;
	}

	void gsl_vector_swap_elements (gsl_vector v, final int i, final int j) {
		double [] data = v.data;

		final int stride = v.stride;
		final int start = v.start;

		if (i != j) {
			double tmp = data[j * stride + start];
			data[j * stride + start] = data[i * stride + start];
			data[i * stride + start] = tmp;

		}
	}

	gsl_vector_view gsl_vector_subvector (gsl_vector v, final int offset, int n) throws Exception {

		if (n == 0) {
			throw new Exception("vector length n must be positive integer");
		}

		if (offset + (n - 1) >= v.size) {
			throw new Exception("view would extend past end of vector");
		}

		gsl_vector s = new gsl_vector();
		gsl_vector_view view = new gsl_vector_view();

		s.data = v.data; // + v.stride * offset;
		s.start = v.start + v.stride * offset;
		s.size = n;
		s.stride = v.stride;

		view.vector = s;
		return view;
	}

	enum CBLAS_UPLO {CblasUpper, CblasLower}

	enum CBLAS_ORDER {CblasRowMajor, CblasColMajor}

	int OFFSET(int N, int incX) {
		return ((incX) > 0 ?  0 : ((N) - 1) * (-(incX)));
	}

	// DNRM2 returns the euclidean norm of a vector
	double cblas_dnrm2 (final int N, final double[] X, final int incX, final int startX) {
		double scale = 0.0;
		double ssq = 1.0;
		int i;
		int ix = startX;

		if (N <= 0 || incX <= 0) {
			return 0;
		} else if (N == 1) {
			return Math.abs(X[0 + startX]);
		}

		for (i = 0; i < N; i++) {
			final double x = X[ix];

			if (x != 0.0) {
				final double ax = Math.abs(x);

				if (scale < ax) {
					ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
					scale = ax;
				} else {
					ssq += (ax / scale) * (ax / scale);
				}
			}

			ix += incX;
		}

		return scale * Math.sqrt(ssq);
	}

	void cblas_dscal (final int N, final double alpha, double [] X, final int incX, final int start) {
		int i;
		int ix;

		if (incX <= 0) {
			return;
		}

		ix = OFFSET(N, incX);

		for (i = 0; i < N; i++) {
			X[start + ix] *= alpha;
			ix += incX;
		}

	}

	void cblas_dsyr2 (final CBLAS_ORDER order, final CBLAS_UPLO Uplo, final int N, final double alpha, final double [] X,
			final int incX, final int startX, final double [] Y, final int incY, final int startY, double [] A, final int lda, final int startA) throws Exception {

		int i, j;

		if (N == 0)
			return;

		if (alpha == 0.0)
			return;

		//		if ((order == CBLAS_ORDER.CblasRowMajor && Uplo == CBLAS_UPLO.CblasUpper)
		//				|| (order == CBLAS_ORDER.CblasColMajor && Uplo == CBLAS_UPLO.CblasLower)) {
		//			int ix = OFFSET(N, incX) + startX;
		//			int iy = OFFSET(N, incY) + startY;
		//			for (i = 0; i < N; i++) {
		//				final double tmp1 = alpha * X[ix];
		//				final double tmp2 = alpha * Y[iy];
		//				int jx = ix;
		//				int jy = iy;
		//				for (j = i; j < N; j++) {
		//					A[lda * i + j + startA] += tmp1 * Y[jy] + tmp2 * X[jx];
		//					jx += incX;
		//					jy += incY;
		//				}
		//				ix += incX;
		//				iy += incY;
		//			}
		//		} else
		if ((order == CBLAS_ORDER.CblasRowMajor && Uplo == CBLAS_UPLO.CblasLower)
				|| (order == CBLAS_ORDER.CblasColMajor && Uplo == CBLAS_UPLO.CblasUpper)) {
			int ix = OFFSET(N, incX) + startX;
			int iy = OFFSET(N, incY) + startY;
			for (i = 0; i < N; i++) {
				final double tmp1 = alpha * X[ix];
				final double tmp2 = alpha * Y[iy];
				int jx = OFFSET(N, incX) + startX;
				int jy = OFFSET(N, incY) + startY;
				for (j = 0; j <= i; j++) {
					A[lda * i + j + startA] += tmp1 * Y[jy] + tmp2 * X[jx];
					jx += incX;
					jy += incY;
				}
				ix += incX;
				iy += incY;
			}
		} else {
			throw new Exception("unrecognized operation");
		}

	}

	void cblas_dsymv (final CBLAS_ORDER order, final CBLAS_UPLO Uplo, final int N, final double alpha, final double[] A, final int lda, final int startA,
			final double[] X, final int incX, final int startX, final double beta, double [] Y, final int incY, final int startY) throws Exception {

		int i, j;

		if (alpha == 0.0 && beta == 1.0)
			return;

		/* form  y := beta*y */
		if (beta == 0.0) {
			int iy = OFFSET(N, incY);
			for (i = 0; i < N; i++) {
				Y[iy + startY] = 0.0;
				iy += incY;
			}
		} else if (beta != 1.0) {
			int iy = OFFSET(N, incY);
			for (i = 0; i < N; i++) {
				Y[iy + startY] *= beta;
				iy += incY;
			}
		}

		if (alpha == 0.0)
			return;

		/* form  y := alpha*A*x + y */

		//		if ((order == CBLAS_ORDER.CblasRowMajor && Uplo == CBLAS_UPLO.CblasUpper)
		//				|| (order == CBLAS_ORDER.CblasColMajor && Uplo == CBLAS_UPLO.CblasLower)) {
		//			int ix = OFFSET(N, incX) + startX;
		//			int iy = OFFSET(N, incY) + startY;
		//			for (i = 0; i < N; i++) {
		//				double temp1 = alpha * X[ix];
		//				double temp2 = 0.0;
		//				final int j_min = i + 1;
		//				final int j_max = N;
		//				int jx = OFFSET(N, incX) + j_min * incX + startX;
		//				int jy = OFFSET(N, incY) + j_min * incY + startY;
		//				Y[iy] += temp1 * A[lda * i + i + startA];
		//				for (j = j_min; j < j_max; j++) {
		//					Y[jy] += temp1 * A[lda * i + j + startA];
		//					temp2 += X[jx] * A[lda * i + j + startA];
		//					jx += incX;
		//					jy += incY;
		//				}
		//				Y[iy] += alpha * temp2;
		//				ix += incX;
		//				iy += incY;
		//			}
		//		} else
		if ((order == CBLAS_ORDER.CblasRowMajor && Uplo == CBLAS_UPLO.CblasLower)
				|| (order == CBLAS_ORDER.CblasColMajor && Uplo == CBLAS_UPLO.CblasUpper)) {
			int ix = OFFSET(N, incX) + (N - 1) * incX + startX;
			int iy = OFFSET(N, incY) + (N - 1) * incY + startY;
			for (i = N; i > 0 && (i-- != 0);) {
				double temp1 = alpha * X[ix];
				double temp2 = 0.0;
				final int j_min = 0;
				final int j_max = i;
				int jx = OFFSET(N, incX) + j_min * incX + startX;
				int jy = OFFSET(N, incY) + j_min * incY + startY;
				Y[iy] += temp1 * A[lda * i + i + startA];
				for (j = j_min; j < j_max; j++) {
					Y[jy] += temp1 * A[lda * i + j + startA];
					temp2 += X[jx] * A[lda * i + j+ startA];
					jx += incX;
					jy += incY;
				}
				Y[iy] += alpha * temp2;
				ix -= incX;
				iy -= incY;
			}
		} else {
			throw new Exception("unrecognized operation");
		}
	}

	double gsl_blas_dnrm2 (final gsl_vector X) {
		return cblas_dnrm2 (X.size, X.data, X.stride, X.start);
	}

	void gsl_blas_dsymv (CBLAS_UPLO Uplo, double alpha, final gsl_matrix A, final gsl_vector X, double beta, gsl_vector Y) throws Exception {
		final int M = A.size1;
		final int N = A.size2;

		if (M != N) {
			throw new Exception("matrix must be square");
		} else if (N != X.size || N != Y.size) {
			throw new Exception("invalid length");
		}

		cblas_dsymv (CBLAS_ORDER.CblasRowMajor, Uplo, N, alpha, A.data, A.tda, A.start, X.data, X.stride, X.start, beta, Y.data, Y.stride, Y.start);
	}

	void gsl_blas_dsyr2 (CBLAS_UPLO Uplo, double alpha, final gsl_vector X, final gsl_vector  Y, gsl_matrix  A) throws Exception {
		final int M = A.size1;
		final int N = A.size2;

		if (M != N) {
			throw new Exception("matrix must be square");
		} else if (X.size != N || Y.size != N) {
			throw new Exception("invalid length");
		}

		cblas_dsyr2 (CBLAS_ORDER.CblasRowMajor, Uplo, N, alpha, X.data, X.stride, X.start, Y.data, Y.stride, Y.start, A.data, A.tda, A.start);
	}

	void gsl_blas_dscal (double alpha, gsl_vector X) {
		cblas_dscal(X.size, alpha, X.data, X.stride, X.start);
	}

	void gsl_blas_daxpy (double alpha, final gsl_vector X, gsl_vector Y) throws Exception {
		if (X.size != Y.size) {
			throw new Exception("invalid length");
		}

		final int N = X.size;
		final double [] Xd = X.data;
		final int incX = X.stride;
		double [] Yd = Y.data;
		final int incY = Y.stride;
		final int startX = X.start;
		final int startY = Y.start;

		int i;

		if (alpha == 0.0) {
			return;
		}

		int ix = OFFSET(N, incX) + startX;
		int iy = OFFSET(N, incY) + startY;

		for (i = 0; i < N; i++) {
			Yd[iy] += alpha * Xd[ix];
			ix += incX;
			iy += incY;
		}

	}

	double gsl_blas_ddot (final gsl_vector X, final gsl_vector Y) throws Exception {
		if (X.size != Y.size) {
			throw new Exception("invalid length");
		}

		final int N = X.size;
		final double [] Xd = X.data;
		final int incX = X.stride;
		final int startX = X.start;
		final double [] Yd = Y.data;
		final int incY = Y.stride;
		final int startY = Y.start;

		double r = 0.0;
		int i;

		int ix = OFFSET(N, incX) + startX;
		int iy = OFFSET(N, incY) + startY;

		for (i = 0; i < N; i++) {
			r += Xd[ix] * Yd[iy];
			ix += incX;
			iy += incY;
		}

		return r;
	}

	double gsl_linalg_householder_transform (gsl_vector v) throws Exception {
		/* replace v[0:n-1] with a householder vector (v[0:n-1]) and
	     coefficient tau that annihilate v[1:n-1] */

		final int n = v.size ;

		if (n == 1)
		{
			return 0.0; /* tau = 0 */
		}
		else
		{ 
			double alpha, beta, tau ;

			gsl_vector_view x = gsl_vector_subvector (v, 1, n - 1) ; 

			double xnorm = gsl_blas_dnrm2 (x.vector);

			if (xnorm == 0) 
			{
				return 0.0; /* tau = 0 */
			}

			alpha = gsl_vector_get (v, 0) ;
			beta = - (alpha >= 0.0 ? +1.0 : -1.0) * Math.hypot(alpha, xnorm) ;
			tau = (beta - alpha) / beta ;

			gsl_blas_dscal (1.0 / (alpha - beta), x.vector);
			gsl_vector_set (v, 0, beta) ;

			return tau;
		}
	}

	void gsl_linalg_householder_hm (double tau, final gsl_vector v, gsl_matrix A) throws Exception {
		/* applies a householder transformation v,tau to matrix m */

		if (tau == 0.0) {
			return;
		}

		/* blas */
		if (true) {
		gsl_vector_view v1 = gsl_vector_subvector (v, 1, v.size - 1);
		gsl_matrix_view A1 = gsl_matrix_submatrix (A, 1, 0, A.size1 - 1, A.size2);
		int j;

		for (j = 0; j < A.size2; j++)
		{
			double wj = 0.0;
			gsl_vector_view A1j = gsl_matrix_column(A1.matrix, j);
			wj = gsl_blas_ddot (A1j.vector, v1.vector);
			wj += gsl_matrix_get(A,0,j);

			{
				double A0j = gsl_matrix_get (A, 0, j);
				gsl_matrix_set (A, 0, j, A0j - tau *  wj);
			}

			gsl_blas_daxpy (-tau * wj, v1.vector, A1j.vector);
		}
		} else { /* no blas */
			int i, j;

			for (j = 0; j < A.size2; j++) {
				/* Compute wj = Akj vk */

				double wj = gsl_matrix_get(A,0,j);  

				for (i = 1; i < A.size1; i++) { /* note, computed for v(0) = 1 above */
					wj += gsl_matrix_get(A,i,j) * gsl_vector_get(v,i);
				}

				/* Aij = Aij - tau vi wj */

				/* i = 0 */
				{
					double A0j = gsl_matrix_get (A, 0, j);
					gsl_matrix_set (A, 0, j, A0j - tau *  wj);
				}

				/* i = 1 .. M-1 */

				for (i = 1; i < A.size1; i++) {
					double Aij = gsl_matrix_get (A, i, j);
					double vi = gsl_vector_get (v, i);
					gsl_matrix_set (A, i, j, Aij - tau * vi * wj);
				}
			}
		}
	}

	/* Factorise a symmetric matrix A into
	 *
	 * A = Q T Q'
	 *
	 * where Q is orthogonal and T is symmetric tridiagonal.  Only the
	 * diagonal and lower triangular part of A is referenced and modified.
	 *
	 * On exit, T is stored in the diagonal and first subdiagonal of
	 * A. Since T is symmetric the upper diagonal is not stored.
	 *
	 * Q is stored as a packed set of Householder transformations in the
	 * lower triangular part of the input matrix below the first subdiagonal.
	 *
	 * The full matrix for Q can be obtained as the product
	 *
	 *       Q = Q_1 Q_2 ... Q_(N-2)
	 *
	 * where 
	 *
	 *       Q_i = (I - tau_i * v_i * v_i')
	 *
	 * and where v_i is a Householder vector
	 *
	 *       v_i = [0, ... , 0, 1, A(i+1,i), A(i+2,i), ... , A(N,i)]
	 *
	 * This storage scheme is the same as in LAPACK.  See LAPACK's
	 * ssytd2.f for details.
	 *
	 * See Golub & Van Loan, "Matrix Computations" (3rd ed), Section 8.3 
	 *
	 * Note: this description uses 1-based indices. The code below uses
	 * 0-based indices 
	 */
	
	void gsl_linalg_symmtd_decomp (gsl_matrix A, gsl_vector tau) throws Exception {
		if (A.size1 != A.size2) {
			throw new Exception("symmetric tridiagonal decomposition requires square matrix");
		} else if (tau.size + 1 != A.size1) {
			throw new Exception("size of tau must be (matrix size - 1)");
		} else {
			final int N = A.size1;
			int i;

			for (i = 0 ; i < N - 2; i++)
			{
				gsl_vector_view c = gsl_matrix_column (A, i);
				gsl_vector_view v = gsl_vector_subvector (c.vector, i + 1, N - (i + 1));
				double tau_i = gsl_linalg_householder_transform (v.vector);

				/* Apply the transformation H^T A H to the remaining columns */

				if (tau_i != 0.0) 
				{
					gsl_matrix_view m = gsl_matrix_submatrix (A, i + 1, i + 1, 
							N - (i+1), N - (i+1));
					double ei = gsl_vector_get(v.vector, 0);
					gsl_vector_view x = gsl_vector_subvector (tau, i, N-(i+1));
					gsl_vector_set (v.vector, 0, 1.0);

					/* x = tau * A * v */
					gsl_blas_dsymv (CBLAS_UPLO.CblasLower, tau_i, m.matrix, v.vector, 0.0, x.vector);

					/* w = x - (1/2) tau * (x' * v) * v  */
					{
						double xv, alpha;
						xv = gsl_blas_ddot(x.vector, v.vector/*, xv*/);
						alpha = - (tau_i / 2.0) * xv;
						gsl_blas_daxpy(alpha, v.vector, x.vector);
					}

					/* apply the transformation A = A - v w' - w v' */
					gsl_blas_dsyr2(CBLAS_UPLO.CblasLower, -1.0, v.vector, x.vector, m.matrix);

					gsl_vector_set (v.vector, 0, ei);
				}

				gsl_vector_set (tau, i, tau_i);
			}

		}
	} 
	
	/*  Form the orthogonal matrix Q from the packed QR matrix */
	void gsl_linalg_symmtd_unpack (final gsl_matrix A, final gsl_vector tau, gsl_matrix Q, gsl_vector diag,	gsl_vector sdiag) throws Exception {
		if (A.size1 !=  A.size2) {
			throw new Exception("matrix A must be square");
		} else if (tau.size + 1 != A.size1) {
			throw new Exception("size of tau must be (matrix size - 1)");
		} else if (Q.size1 != A.size1 || Q.size2 != A.size1) {
			throw new Exception("size of Q must match size of A");
		} else if (diag.size != A.size1) {
			throw new Exception("size of diagonal must match size of A");
		} else if (sdiag.size + 1 != A.size1) {
			throw new Exception("size of subdiagonal must be (matrix size - 1)");
		} else {
			final int N = A.size1;

			int i;

			/* Initialize Q to the identity */

			gsl_matrix_set_identity (Q);

			for (i = N - 2; i-- > 0;)
			{
				gsl_vector_view c = gsl_matrix_column (A, i);
				gsl_vector_view h = gsl_vector_subvector (c.vector, i + 1, N - (i+1));
				double ti = gsl_vector_get (tau, i);

				gsl_matrix_view m = gsl_matrix_submatrix (Q, i + 1, i + 1, N-(i+1), N-(i+1));

				gsl_linalg_householder_hm (ti, h.vector, m.matrix);
			}

			/* Copy diagonal into diag */

			for (i = 0; i < N; i++)
			{
				double Aii = gsl_matrix_get (A, i, i);
				gsl_vector_set (diag, i, Aii);
			}

			/* Copy subdiagonal into sd */

			for (i = 0; i < N - 1; i++)
			{
				double Aji = gsl_matrix_get (A, i+1, i);
				gsl_vector_set (sdiag, i, Aji);
			}

		}
	}


	/* remove off-diagonal elements which are neglegible compared with the
	   neighboring diagonal elements */
	void chop_small_elements (final int N, final double d[], double sd[]) {
		double d_i = d[0];

		int i;

		for (i = 0; i < N - 1; i++)
		{
			double sd_i = sd[i];
			double d_ip1 = d[i + 1];

			double GSL_DBL_EPSILON = 2.2204460492503131e-16;

			if (Math.abs (sd_i) < GSL_DBL_EPSILON * (Math.abs (d_i) + Math.abs (d_ip1)))
			{
				sd[i] = 0.0;
			}
			d_i = d_ip1;
		}
	}

	/* Generate a Givens rotation (cos,sin) which takes v=(x,y) to (|v|,0) 
	   From Golub and Van Loan, "Matrix Computations", Section 5.1.8 */
	double [] /* c, s */ create_givens (final double a, final double b) {
		double c, s;

		if (b == 0)
		{
			c = 1;
			s = 0;
		}
		else if (Math.abs (b) > Math.abs (a))
		{
			double t = -a / b;
			double s1 = 1.0 / Math.sqrt (1 + t * t);
			s = s1;
			c = s1 * t;
		}
		else
		{
			double t = -b / a;
			double c1 = 1.0 / Math.sqrt (1 + t * t);
			c = c1;
			s = c1 * t;
		}

		double[] r = {c,s};

		return r;
	}

	double trailing_eigenvalue (final int n, final double d[], final double sd[], int a) {
		double ta = d[n - 2 + a];
		double tb = d[n - 1 + a];
		double tab = sd[n - 2 + a];

		double dt = (ta - tb) / 2.0;

		double mu;

		if (dt > 0)
		{
			mu = tb - tab * (tab / (dt + Math.hypot (dt, tab)));
		}
		else if (dt == 0) 
		{
			mu = tb - Math.abs(tab);
		}
		else
		{
			mu = tb + tab * (tab / ((-dt) + Math.hypot (dt, tab)));
		}

		return mu;
	}

	void qrstep (final int n, double d[], double sd[], double gc[], double gs[], int a) {
		double x, z;
		double ak, bk, zk, ap, bp, aq, bq;
		int k;

		double mu = trailing_eigenvalue (n, d, sd, a);

		x = d[0 + a] - mu;
		z = sd[0 + a];

		ak = 0;
		bk = 0;
		zk = 0;

		ap = d[0 + a];
		bp = sd[0 + a];

		aq = d[1 + a];

		if (n == 2)
		{

			double[] r = create_givens (x, z);
			double c = r[0], s = r[1];

			if (gc != null)
				gc[0] = c; 
			if (gs != null)
				gs[0] = s;

			{
				double ap1 = c * (c * ap - s * bp) + s * (s * aq - c * bp);
				double bp1 = c * (s * ap + c * bp) - s * (s * bp + c * aq);

				double aq1 = s * (s * ap + c * bp) + c * (s * bp + c * aq);

				ak = ap1;
				bk = bp1;

				ap = aq1;
			}

			d[0 + a] = ak;
			sd[0 + a] = bk;
			d[1 + a] = ap;

			return;
		}

		bq = sd[1 + a];

		for (k = 0; k < n - 1; k++)
		{
			double[] r = create_givens (x, z);
			double c = r[0], s = r[1];

			/* store Givens rotation */
			if (gc != null)
				gc[k] = c; 
			if (gs != null)
				gs[k] = s;

			/* compute G' T G */

			{
				double bk1 = c * bk - s * zk;

				double ap1 = c * (c * ap - s * bp) + s * (s * aq - c * bp);
				double bp1 = c * (s * ap + c * bp) - s * (s * bp + c * aq);
				double zp1 = -s * bq;

				double aq1 = s * (s * ap + c * bp) + c * (s * bp + c * aq);
				double bq1 = c * bq;

				ak = ap1;
				bk = bp1;
				zk = zp1;

				ap = aq1;
				bp = bq1;

				if (k < n - 2)
					aq = d[k + 2 + a];
				if (k < n - 3)
					bq = sd[k + 2 + a];

				d[k + a] = ak;

				if (k > 0)
					sd[k - 1 + a] = bk1;

				if (k < n - 2)
					sd[k + 1 + a] = bp;

				x = bk;
				z = zk;
			}
		}

		/* k = n - 1 */
		d[k + a] = ap;
		sd[k - 1 + a] = bk;
	}

	/* Compute eigenvalues/eigenvectors of real symmetric matrix using
	   reduction to tridiagonal form, followed by QR iteration with
	   implicit shifts.

	   See Golub & Van Loan, "Matrix Computations" (3rd ed), Section 8.3
	   */
	
	void gsl_eigen_symmv (gsl_matrix A, gsl_vector eval, gsl_matrix evec, gsl_eigen_symmv_workspace w) throws Exception {
		if (A.size1 != A.size2) {
			throw new Exception("matrix must be square to compute eigenvalues");
		} else if (eval.size != A.size1) {
			throw new Exception("eigenvalue vector must match matrix size");
		} else if (evec.size1 != A.size1 || evec.size2 != A.size1) {
			throw new Exception("eigenvector matrix must match matrix size");
		} else {
			final double [] d = w.d;
			final double [] sd = w.sd;
			final int N = A.size1;
			int a, b;

			/* handle special case */

			if (N == 1)
			{
				double A00 = gsl_matrix_get (A, 0, 0);
				gsl_vector_set (eval, 0, A00);
				gsl_matrix_set (evec, 0, 0, 1.0);
				return;
			}

			/* use sd as the temporary workspace for the decomposition when
	         computing eigenvectors */

			{
				gsl_vector_view d_vec = new gsl_vector_view (d, N);
				gsl_vector_view sd_vec = new gsl_vector_view (sd, N - 1);
				gsl_vector_view tau = new gsl_vector_view (sd, N - 1);
				gsl_linalg_symmtd_decomp (A, tau.vector);
				gsl_linalg_symmtd_unpack (A, tau.vector, evec, d_vec.vector, sd_vec.vector);
			}

			/* Make an initial pass through the tridiagonal decomposition
	         to remove off-diagonal elements which are effectively zero */

			chop_small_elements (N, d, sd);

			/* Progressively reduce the matrix until it is diagonal */

			b = N - 1;

			while (b > 0)
			{
				if (sd[b - 1] == 0.0 || Double.isNaN(sd[b - 1]))
				{
					b--;
					continue;
				}

				/* Find the largest unreduced block (a,b) starting from b
	             and working backwards */

				a = b - 1;

				while (a > 0)
				{
					if (sd[a - 1] == 0.0)
					{
						break;
					}
					a--;
				}

				{
					int i;
					final int n_block = b - a + 1;
					double[] gc = w.gc;
					double[] gs = w.gs;

					/* apply QR reduction with implicit deflation to the unreduced block */

					qrstep (n_block, d, sd, gc, gs, a);

					/* Apply  Givens rotation Gij(c,s) to matrix Q,  Q <- Q G */

					for (i = 0; i < n_block - 1; i++)
					{
						final double c = gc[i], s = gs[i];
						int k;

						for (k = 0; k < N; k++)
						{
							double qki = gsl_matrix_get (evec, k, a + i);
							double qkj = gsl_matrix_get (evec, k, a + i + 1);
							gsl_matrix_set (evec, k, a + i, qki * c - qkj * s);
							gsl_matrix_set (evec, k, a + i + 1, qki * s + qkj * c);
						}
					}

					/* remove any small off-diagonal elements */

					chop_small_elements (N, d, sd);
				}
			}

			{
				gsl_vector_view d_vec = new gsl_vector_view (d, N);

				for (int i = 0; i < N; i++) {
					gsl_vector_set(eval, i, gsl_vector_get(d_vec.vector, i));
				}
			}

			return;
		}
	}


	void gsl_eigen_symmv_sort (gsl_vector eval, gsl_matrix evec) throws Exception {

		if (evec.size1 != evec.size2) {
			throw new Exception("eigenvector matrix must be square");
		} else if (eval.size != evec.size1) {
			throw new Exception("eigenvalues must match eigenvector matrix");
		} else {
			final int N = eval.size;
			int i;

			for (i = 0; i < N - 1; i++)
			{
				int j;
				int k = i;

				double ek = gsl_vector_get (eval, i);

				/* search for something to swap */
				for (j = i + 1; j < N; j++) {
					boolean test;
					final double ej = gsl_vector_get (eval, j);

					/*      switch (sort_type)
	                {       
	                case GSL_EIGEN_SORT_VAL_ASC:
	                  test = (ej < ek);
	                  break;
	                case GSL_EIGEN_SORT_VAL_DESC:
	                  test = (ej > ek);
	                  break;
	                case GSL_EIGEN_SORT_ABS_ASC:
	                  test = (fabs (ej) < fabs (ek));
	                  break;
	                case GSL_EIGEN_SORT_ABS_DESC:*/
					test = (Math.abs (ej) > Math.abs (ek));
					/*   break;
	                default:
	                  throw new Exception("unrecognized sort type");
	                }*/

					if (test)
					{
						k = j;
						ek = ej;
					}
				}

				if (k != i)
				{
					/* swap eigenvalues */
					gsl_vector_swap_elements (eval, i, k);

					/* swap eigenvectors */
					gsl_matrix_swap_columns (evec, i, k);
				}
			}
		}
	}
}
