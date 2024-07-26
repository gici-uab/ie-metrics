package GiciMatrix;

/**
 * A set of tools to operate on matrices.
 *
 * <p> Matrices are represented using an array of arrays of floats.
 * There are various kinds of matrices: complete, upper Hessenberg, upper triangular, etc.
 * The difference between all of them is the row length, which may be variable.
 *
 * <p> For example a triangular matrix may have the following structure:
 * <blockquote>
 * float[][] matrix = new float[2][];
 * matrix[0] = new float[2];
 * matrix[1] = new float[1];
 * </blockquote>
 *
 * <p> Methods that only work with one kind of matrices have a name that ends in XXYY where
 * XX and YY are either C (complete), UH (upper Hessenberg) or UT (upper triangular). If a method
 * has one input it will end in XX, if it has two then XXYY, and so on.
 *
 * @author Ian Blanes
 *
 */
public class MatrixAlgebra {
	/**
	 * This is an utility class and shall not be constructed.
	 */
	protected MatrixAlgebra() {
		throw new UnsupportedOperationException();
	}

	/**
	 * Transpose of a Matrix. Only works on complete matrices.
	 *
	 * @param a the complete matrix to be transposed.
	 * @return the transpose of a.
	 */
	public static float[][] transposeC(final float[][] a) {
		float[][] r = new float[a[0].length][a.length];

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[i].length; j++) {
				assert (a[i].length == a[0].length);
				r[j][i] = a[i][j];
			}
		}

		return r;
	}

	/**
	 * Multiplication extending one matrix with the identity.
	 *
	 * @param a is extended with the identity on the upper left side and zeros
	 * elsewhere to match the size of b. A must be complete and square.
	 * @param b a complete matrix.
	 * @return the result.
	 */
	public static float[][] subMultiplicationCC(final float [][] a, final float[][] b) {
		assert (a.length > 0 && a.length == a[0].length);
		assert (a.length <= b.length);
		assert (b.length > 0);

		float[][] r = new float[b.length][b[0].length];

		int skip = b.length - a.length;
		int i, j, k;

		/* Calculate the zone the height of a */
		for (i = 0; i < a.length; i++) {
			assert (a.length == a[i].length);
			assert (b[0].length == b[i + skip].length);

			for (j = 0; j < b[0].length; j++) {
				for (k = 0; k < a.length; k++) {
					r[i + skip][j] += a[i][k] * b[k + skip][j];
				}
			}
		}

		/* Fill the rest with b */
		for (i = 0; i < skip; i++) {
			assert (b[0].length == b[i].length);

			for (j = 0; j < b[0].length; j++) {
				r[i][j] = b[i][j];
			}
		}

		return r;
	}

	/**
	 * Complete matrix multiplication.
	 *
	 * @param a complete matrix.
	 * @param b complete matrix.
	 * @return A x B.
	 */
	public static float[][] multiplicationCC(final float [][] a, final float[][] b) {
		assert (a.length > 0);
		assert (b.length > 0);
		assert (a[0].length == b.length);

		float[][] r = new float[a.length][b[0].length];

		for (int i = 0; i < a.length; i++) {
			assert (a[0].length == a[i].length);

			for (int j = 0; j < b[0].length; j++) {
				double acc = 0.0;
				
				for (int k = 0; k < b.length; k++) {
					assert (b[0].length == b[k].length);

					acc += a[i][k] * b[k][j];
				}
				
				r[i][j] = (float) acc;
			}
		}

		return r;
	}

	/**
	 * Matrix multiplication between a upper triangular matrix and a upper Hessenberg matrix.
	 *
	 * @param a upper triangular matrix.
	 * @param b upper Hessenberg matrix.
	 * @return A x B.
	 */
	public static float[][] multiplicationUTUH(final float [][] a, final float[][] b) {

		assert (a.length > 0 && a.length == a[0].length);
		assert (b.length > 0 && b.length == b[0].length);
		assert (a.length == b.length);

		float[][] r = zeroUH(a.length);

		for (int i = 0; i < r.length; i++) {
			int skipR = Math.max(0, i - 1);

			for (int j = skipR; j < r.length; j++) {
				int skipA = i;

				for (int k = skipA; k < Math.min(j + 2, b.length); k++) {
					int skipB = Math.max(0, k - 1);

					r[i][j - skipR] += a[i][k - skipA] * b[k][j - skipB];
				}
			}
		}

		return r;
	}

	/**
	 * Matrix multiplication between a complete matrix and a upper Hessenberg matrix.
	 *
	 * @param a complete matrix.
	 * @param b upper Hessenberg matrix.
	 * @return A x B.
	 */
	public static float[][] multiplicationCUH(final float [][] a, final float[][] b) {
		assert (b.length > 0 && b.length == b[0].length);

		/* Check for compatibility */
		assert (a.length > 0 && a[0].length == b.length);

		float[][] r = new float[a.length][b.length];

		for (int i = 0; i < a.length; i++) {
			assert (a[0].length == a[i].length);

			for (int j = 0; j < b.length; j++) {
				int uhLength = Math.min(j + 2, b.length);

				/* Check for UH */
				assert (b[j].length == b.length - Math.max(j - 1, 0));

				for (int k = 0; k < uhLength; k++) {
					int skip = Math.max(0, k - 1);
					r[i][j] += a[i][k] * b[k][j - skip];
				}
			}
		}

		return r;
	}

	/**
	 * Matrix-Vector multiplication between a complete matrix and one vector.
	 *
	 * @param a complete matrix.
	 * @param v column vector.
	 * @return A x v.
	 */
	public static float[] multiplicationCV(final float [][] a, final float [] v) {
		assert (a.length > 0 && a[0].length == v.length);

		float[] r = new float[a.length];

		for (int i = 0; i < a.length; i++) {
			assert (a[0].length == a[i].length);

			for (int j = 0; j < v.length; j++) {
				r[i] += a[i][j] * v[j];
			}
		}

		return r;
	}

	/**
	 * Complete identity matrix.
	 *
	 * @param n size.
	 * @return Id(n).
	 */
	public static float[][] identityC(final int n) {
		assert (n > 0);

		float [][] r = new float[n][n];

		for (int i = 0; i < n; i++) {
			r[i][i] = 1.0f;
		}

		return r;
	}

	/**
	 * Upper triangular identity matrix.
	 *
	 * @param n size.
	 * @return Id(n).
	 */
	public static float[][] identityUT(final int n) {
		assert (n > 0);

		float [][] r = new float[n][];

		for (int i = 0; i < n; i++) {
			r[i] = new float[n - i];
			r[i][0] = 1.0f;
		}

		return r;
	}

	/**
	 * Upper Hessenberg matrix filled with zeros.
	 *
	 * @param n size.
	 * @return Zero(n).
	 */
	public static float[][] zeroUH(final int n) {
		assert (n > 0);

		float [][] r = new float[n][];

		r[0] = new float[n];

		for (int i = 1; i < n; i++) {
			r[i] = new float[n - i + 1];
		}

		return r;
	}

	/**
	 * Upper Hessenberg identity matrix.
	 *
	 * @param n size.
	 * @return Id(n).
	 */
	public static float[][] identityUH(final int n) {
		assert (n > 0);

		float [][] r = zeroUH(n);

		r[0][0] = 1.0f;

		for (int i = 1; i < n; i++) {
			r[i][1] = 1.0f;
		}

		return r;
	}

	/**
	 * Test for equality. Assumes that both matrices are aligned properly
	 * and have the same sparsity.
	 *
	 * @param a a matrix.
	 * @param b another matrix.
	 * @param epsilon the maximum difference that will be allowed and wont cause matrices to differ.
	 * @return a boolean stating whether the matrices are equal(<code>true</code>) or not(<code>false</code>).
	 */
	public static boolean compare(final float[][] a, final float[][] b, final float epsilon) {

		assert (a.length == b.length);

		for (int i = 0; i < a.length; i++) {
			assert (a[i].length == b[i].length);

			for (int j = 0; j < a[i].length; j++) {
				if (Math.abs(a[i][j] - b[i][j]) > epsilon) {
					return false;
				}
			}
		}

		return true;
	}

	/**
	 * Cuts an upper Hessenberg matrix into a upper triangular. Elements outside the
	 * the resulting matrix are deprecated.
	 *
	 * @param a an upper Hessenberg matrix.
	 * @return an upper triangular matrix.
	 */
	public static float[][] cutUHUT(final float[][] a) {
		float [][] r = identityUT(a.length);

		for (int i = 0; i < a.length; i++) {
			/* Test for a upper Hessenberg */
			assert (i == 0 ? a[i].length == a.length : a[i].length + i - 1 == a.length);

			r[i] = new float[a.length - i];
			System.arraycopy(a[i], (i == 0 ? 0 : 1), r[i], 0, a.length - i);
		}

		return r;
	}

	/**
	 * Cuts a complete matrix into a upper Hessenberg. Elements outside the
	 * the resulting matrix are deprecated.
	 *
	 * @param a a complete matrix.
	 * @return an upper Hessenberg matrix.
	 */
	public static float[][] cutCUH(final float[][] a) {
		assert (a.length > 0 && a[0].length == a.length);

		float [][] r = identityUH(a.length);
		r[0] = new float[a.length];
		System.arraycopy(a[0], 0, r[0], 0, a.length);

		for (int i = 1; i < a.length; i++) {
			/* Test for a complete matrix */
			assert (a[i].length == a.length);

			r[i] = new float[a.length - i + 1];
			System.arraycopy(a[i], i - 1, r[i], 0, r[i].length);
		}

		return r;
	}

	/**
	 * Copies a matrix.
	 *
	 * @param a input matrix.
	 * @return a copy of a.
	 */
	public static float[][] copy(final float[][] a) {
		float [][] r = new float[a.length][];

		for (int i = 0; i < a.length; i++) {
			r[i] = new float[a[i].length];
			System.arraycopy(a[i], 0, r[i], 0, a[i].length);
		}

		return r;
	}

	/**
	 * Computes the difference between b and a (a - b). Assumes that both matrices are aligned properly
	 * and have the same sparsity.
	 *
	 * @param a a matrix.
	 * @param b another matrix.
	 * @return A - B.
	 */
	public static float[][] substract(final float[][] a, final float[][] b) {
		float[][] r = new float[a.length][];

		assert (a.length == b.length);

		for (int i = 0; i < a.length; i++) {
			assert (a[i].length == b[i].length);

			r[i] = new float[a[i].length];

			for (int j = 0; j < a[i].length; j++) {
				r[i][j] = a[i][j] - b[i][j];
			}
		}

		return r;
	}

	/* DEBUG */
	/**
	 * Fills a matrix with zeros till is a complete matrix.
	 * The width of the resulting matrix is assumed to be the width of the longest row.
	 * Rows are assumed to be aligned on the right, so are padded on the left.
	 *
	 * @param a a matrix.
	 * @return the matrix a padded with zeros.
	 */
	public static float[][] toComplete(final float[][] a) {
		/* Assume it's right aligned */
		int i, width = 0;

		for (i = 0; i < a.length; i++) {
			width = Math.max(width, a[i].length);
		}

		float[][] r = new float[a.length][width];

		for (i = 0; i < a.length; i++) {
			for (int j = 0; j < a[i].length; j++) {
				r[i][j + width - a[i].length] = a[i][j];
			}
		}

		return r;
	}
	
	/**
	 * Pads a square complete matrix with the identity in the bottom right zone.
	 *
	 * @param p input matrix.
	 * @param size padding size.
	 * @return padded p.
	 */
	public static float[][] padIdentityBR(final float[][] p, final int size) {
		float[][] r = MatrixAlgebra.identityC(p.length + size);

		/* copy the matrix */
		for (int j = 0; j < p.length; j++) {
			for (int k = 0; k < p.length; k++) {
				r[j][k] = p[j][k];
			}
		}
		
		return r;
	}
	
	/**
	 * Pads a square complete matrix with the identity in the top left zone.
	 *
	 * @param p input matrix.
	 * @param size padding size.
	 * @return padded p.
	 */
	public static float[][] padIdentityTL(final float[][] p, final int size) {
		float[][] r = MatrixAlgebra.identityC(p.length + size);

		/* copy the rest */
		for (int j = 0; j < p.length; j++) {
			for (int k = 0; k < p.length; k++) {
				r[j + size][k + size] = p[j][k];
			}
		}

		return r;
	}

	/**
	 * Computes the relative difference between b and a (a - b) / b. Assumes that both matrices are aligned properly
	 * and have the same sparsity. b_ij = 0 is not checked.
	 *
	 * @param a Matrix b + error.
	 * @param b Original Matrix.
	 * @return the relative error matrix.
	 */
	public static float[][] relativeError(final float[][] a, final float[][] b) {
		float[][] r = new float[a.length][];

		assert (a.length == b.length);

		for (int i = 0; i < a.length; i++) {
			assert (a[i].length == b[i].length);

			r[i] = new float[a[i].length];

			for (int j = 0; j < a[i].length; j++) {
				r[i][j] = Math.abs((a[i][j] - b[i][j]) / b[i][j]);
			}
		}

		return r;
	}

	public static float determinant(final float[][] a) {
		return determinant(a, 1);
	}
	
	public static float determinant(final float[][] a, final int square) {
		assert (a.length > 0 && a[0].length == a.length);
		
		float[][] b = copy(a);
		int N = a.length;
	
		double ldet = 0;
		double detsign = 1;
		
		/* Gaussian elimination */
		for (int i = 0; i < N; i++) {
			assert(a[i].length == N);
			
			/* Compute the permutation if need */
			if (Math.abs(b[i][i]) < Float.MIN_VALUE) {
				// We are going to swap the rows i and i+k
				int k = 1;
				
				while (i + k < N && Math.abs(b[i+k][i]) < Float.MIN_VALUE) {
					k++;
				}
				
				if (i + k >= N) {
					//Non-singular matrix
					return 0;
				}
				
				float[] t;
				
				t = b[i];
				b[i] = b[i+k];
				b[i+k] = t;
			}
			
			for (int j = i + 1; j < N; j++) {
				b[i][j] /= b[i][i];
			}
			
			ldet += Math.log(Math.abs(b[i][i]));
			detsign *= Math.copySign(1, b[i][i]);
			
			b[i][i] = 1;
			
			for (int j = i + 1; j < N; j++) {
				for (int k = i + 1; k < N; k++) {
					b[j][k] -= b[j][i] * b[i][k];
				}
				
				b[j][i] = 0;
			}
		}
		
		// FIXME! oju amb el signe de les arrels
		return (float)Math.copySign(Math.exp(ldet/(double)square), detsign);
	}
	
	/**
	 * Computes a matrix inverse though Gauss-Jordan elimination.
	 * Numerical Stability has been taken into account
	 * @param a
	 * @return
	 */
	public static float[][] inverseC (final float[][] a) {
		
		assert (a.length > 0 && a.length == a[0].length);
		
		float[][] id = copy(a);
		float[][] a1 = identityC(a.length);
		
		/* Perform gauss-jordan elimination on both i and a1 (only row elementary operation allowed):
		 *	1. Interchange two rows.
		 *	2. Multiply a row with a nonzero number.
		 *	3. Add a row to another one multiplied by a number.
		 */
		
		for (int i = 0; i < id.length; i++) {
			// Select the largest element to be the pivot of the i'th col
			
			float max = Float.NEGATIVE_INFINITY;
			int pos = -1;
			
			for (int j = i; j < id.length; j++) {
				// Select the largest element to be the pivot
				if (max <= Math.abs(id[j][i])) {
					max = Math.abs(id[j][i]);
					pos = j;
				}
			}
			
			assert (pos >= 0);
			
			// Swap rows i <-> pos
			float[] t;

			t = id[i];
			id[i] = id[pos];
			id[pos] = t;

			t = a1[i];
			a1[i] = a1[pos];
			a1[pos] = t;

			// Normalize the pivot
			float divisor = id[i][i];
			//System.out.println(divisor);			
			// Singular matrix
			assert(Math.abs(divisor) > 0.0000001f);
			
			for (int j = i; j < id.length; j++) {
				id[i][j] /= divisor;
			}

			for (int j = 0; j < id.length; j++) {
				a1[i][j] /= divisor;
			}
			
			// Zero above and below the pivot

			for (int j = 0; j < i; j++) {
				float rowfactor = id[j][i]; 
				
				for (int k = i + 1; k < id.length; k++) {
					id[j][k] -= rowfactor * id[i][k]; 
				}
				
				for (int k = 0; k < id.length; k++) {
					a1[j][k] -= rowfactor * a1[i][k]; 
				}
			}

			for (int j = i + 1; j < id.length; j++) {
				float rowfactor = id[j][i]; 
				
				for (int k = i + 1; k < id.length; k++) {
					id[j][k] -= rowfactor * id[i][k]; 
				}
				
				for (int k = 0; k < id.length; k++) {
					a1[j][k] -= rowfactor * a1[i][k];
				}
			}
		}
		
		return a1;
	}
	
	public static float[] matrixDiagonal (final float[][] a) {
		assert(a.length > 0 && a.length == a[0].length);
		
		float[] r = new float[a.length];
		
		for (int i = 0; i < a.length; i++) 
			r[i] = a[i][i];
		
		return r;
	}
	
	
	/**
	 * Prints a matrix in the standard output. Fills a matrix with zeros till is a complete matrix.
	 * The width of the resulting matrix is assumed to be the width of the longest row.
	 * Rows are assumed to be aligned on the right, so are padded on the left.
	 *
	 * @param a a matrix.
	 */
	public static void printMatrix(final float [][] a) {
		/* Assume it's right aligned */
		int i, width = 0;

		for (i = 0; i < a.length; i++) {
			width = Math.max(width, a[i].length);
		}

		for (i = 0; i < a.length; i++) {
			int j;
			/* Padd */
			for (j = 0; j < width - a[i].length; j++) {
				System.out.format("[%10.5g] \t", 0.0f);
			}

			/* Values */
			for (j = 0; j < a[i].length; j++) {
				System.out.format("% 11.5g \t", a[i][j]);
			}

			System.out.println();
		}

		System.out.println();
	}
	
	public static void printMatrix(final long [][] a) {
		/* Assume it's right aligned */
		int i, width = 0;

		for (i = 0; i < a.length; i++) {
			width = Math.max(width, a[i].length);
		}

		for (i = 0; i < a.length; i++) {
			int j;
			/* Padd */
			for (j = 0; j < width - a[i].length; j++) {
				System.out.format("[%d] \t", 0);
			}

			/* Values */
			for (j = 0; j < a[i].length; j++) {
				System.out.format("%d \t", a[i][j]);
			}

			System.out.println();
		}

		System.out.println();
	}
	
	public static void printVector(final float [] a) {
		for (int i = 0; i < a.length; i++) {
			System.out.format("% 11.5g \t", a[i]);
		}

		System.out.println();
	}
	
	public static void printVector(final int [] a) {
		for (int i = 0; i < a.length; i++) {
			System.out.format("% 11d \t", a[i]);
		}

		System.out.println();
	}
}
