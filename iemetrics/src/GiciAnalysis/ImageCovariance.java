package GiciAnalysis;

import java.util.Random;
import GiciMatrix.MatrixAlgebra;

/**
 * Generates a covariance matrix.
 *
 * @author Ian Blanes
 */
public class ImageCovariance {
	/**
	 * This is an utility class and shall not be constructed.
	 */
	protected ImageCovariance() {
		throw new UnsupportedOperationException();
	}

	/* warning sampling = 1 does not guaranted all samples are evaluated, as some samples might be sampled twice! */
	public static float[][] generateCovarianceMatrix(final float[][][] imageSamples, final int dimension, final float samplingMeans, final float samplingCovariance) {
		/* first some contractual verifications */
		final int maxDimensions = 3;
		assert (dimension >= 0 && dimension < maxDimensions);
		assert (imageSamples.length > 0 && imageSamples[0].length > 0 && imageSamples[0][0].length > 0);
		assert (samplingMeans > 0 && samplingMeans <= 1);
		assert (samplingCovariance > 0 && samplingCovariance <= 1);

		/* result */
		float[][] r = null;

		int d0 = imageSamples.length;
		int d1 = imageSamples[0].length;
		int d2 = imageSamples[0][0].length;

		int i, j;

		switch (dimension) {
		case 0:
			r = MatrixAlgebra.identityUT(d0);

			/* decide the sampled locations */
			Random rand = new ParkMillerPRNG(13);

			int meanSamples = (int) (d1 * d2 * samplingMeans);
			int covarianceSamples = (int) (d1 * d2 * samplingCovariance);

			int locations = Math.max(meanSamples, covarianceSamples);

			/* its important for stability that both measures have the same samples */
			meanSamples = locations;
			covarianceSamples = locations;
			
			int[] sampleD1 = new int[locations];
			int[] sampleD2 = new int[locations];

			for (i = 0; i < sampleD1.length; i++) {
				sampleD1[i] = rand.nextInt(d1);
				sampleD2[i] = rand.nextInt(d2);
			}

			/* compute the averages */
			double means[] = new double[d0];

			for (i = 0; i < d0; i++) {
				for (j = 0; j < meanSamples; j++) {
					means[i] += imageSamples[i][sampleD1[j]][sampleD2[j]];
				}

				means[i] /= (double) meanSamples;
			}

			/* For each of the destinations */

			for (i = 0; i < d0; i++) {
				for (j = i; j < d0; j++) {
					/* Calculate between d0=i and d0=j */
					double eXY = 0.0f;

					for (int m = 0; m < covarianceSamples; m++) {				
						//eXY += imageSamples[i][sampleD1[m]][sampleD2[m]] * imageSamples[j][sampleD1[m]][sampleD2[m]];
						eXY += (imageSamples[i][sampleD1[m]][sampleD2[m]] - means[i]) * (imageSamples[j][sampleD1[m]][sampleD2[m]] - means[j]);
					}
					
					// This method is almost mandatory to prevent a catastrophic cancellation
					// There is an alternate method using a recurrence function.
					double cov = (eXY / (double) (covarianceSamples - 1));
						//- means[i] * means[j] * (covarianceSamples / (double)(covarianceSamples - 1));

					assert (cov >= 0 || i != j);
					
					r[i][j - i] = (float)cov;
				}
			}

			break;
		default:
			throw new ArrayIndexOutOfBoundsException();
		}

		return r;
	}

	/**
	 * Generates a covariance matrix. This method generates the full and exact covariance matrix.
	 *
	 * @param imageSamples used to compute the covariance.
	 * @param dimension where the covariance is computed over (z=0, y=1, x=2).
	 * @return an upper triangularar matrix with the result (@see MatrixAlgebra).
	 */
	public static float[][] generateCovarianceMatrix(final float[][][] imageSamples, final int dimension) {
		/* first some contractual verifications */
		final int maxDimensions = 3;
		assert (dimension >= 0 && dimension < maxDimensions);
		assert (imageSamples.length > 0);
		assert (imageSamples[0].length > 0);
		assert (imageSamples[0][0].length > 0);

		/* result */
		float[][] r = null;

		int d0 = imageSamples.length;
		int d1 = imageSamples[0].length;
		int d2 = imageSamples[0][0].length;

		int i, j, k, l, n;

		switch (dimension) {
		case 0:
			r = MatrixAlgebra.identityUT(d0);

			n = d1 * d2;

			assert(n > 0);
			
			/* compute the averages */
			double means[] = new double[d0];			

			for (i = 0; i < d0; i++) {
				for (k = 0; k < d1; k++) {
					for (l = 0; l < d2; l++) {
						means[i] += imageSamples[i][k][l];
					}
				}

				means[i] /= (double) (n);
			}

			/* For each of the destinations */
			for (i = 0; i < d0; i++) {
				for (j = i; j < d0; j++) {
					/* Calculate between d0=i and d0=j */
					double eXY = 0.0f;

					for (k = 0; k < d1; k++) {
						for (l = 0; l < d2; l++) {
							//eXY += imageSamples[i][k][l] * imageSamples[j][k][l];
							eXY += (imageSamples[i][k][l] - means[i]) * (imageSamples[j][k][l] - means[j]);
						}
					}
					//warning: catastrophic cancellation occurs using the other method.
					double cov = (eXY / (double) (n - 1)); // - means[i] * means[j] * (n / (double)(n - 1));

					r[i][j - i] = (float)cov;
				}
			}

			break;
		case 1:
			r = MatrixAlgebra.identityUT(d1);

			/* For each of the destinations */
			n = d0 * d2;

			for (i = 0; i < d1; i++) {
				for (j = i; j < d1; j++) {
					/* Calculate between d0=i and d0=j */
					double eXY = 0.0f;
					double eX = 0.0f;
					double eY = 0.0f;

					for (k = 0; k < d0; k++) {
						for (l = 0; l < d2; l++) {
							eXY += imageSamples[k][i][l] * imageSamples[k][j][l];
							eX += imageSamples[k][i][l];
							eY += imageSamples[k][j][l];
						}
					}
					// debugme first!
					assert(false);
					float cov = (float)((eXY / (double) (n - 1)) - (eX / (double) (n - 1)) * (eY / (double) (n - 1)));

					r[i][j - i] = cov;
				}
			}

			break;
		case 2:
			r = MatrixAlgebra.identityUT(d2);

			/* For each of the destinations */
			n = d0 * d1;

			for (i = 0; i < d2; i++) {
				for (j = i; j < d2; j++) {
					/* Calculate between d0=i and d0=j */
					double eXY = 0.0f;
					double eX = 0.0f;
					double eY = 0.0f;

					for (k = 0; k < d0; k++) {
						for (l = 0; l < d1; l++) {
							eXY += imageSamples[k][l][i] * imageSamples[k][l][j];
							eX += imageSamples[k][l][i];
							eY += imageSamples[k][l][j];
						}
					}
					// debugme first!
					assert(false);
					float cov = (float)((eXY / (double) (n - 1)) - (eX / (double) (n - 1)) * (eY / (double) (n - 1)));

					r[i][j - i] = cov;
				}
			}

			break;
		default:
			throw new ArrayIndexOutOfBoundsException();
		}

		return r;
	}

	/**
	 *  
	 * @param cov an upper triangularar matrix with the covariance matrix (@see MatrixAlgebra).
	 * @return
	 */
	public static float[][] generateCorrelationMatrixFromCovarianceMatrix (final float[][] cov) {
		float[][] corr = MatrixAlgebra.copy(cov);

		float[] sigma = new float[corr.length];

		/* get the correct sigmas */
		for (int i = 0; i < corr.length; i++) {
			sigma[i] = (float) Math.sqrt(corr[i][0]);
			assert(! Float.isNaN(sigma[i]));
		}

		/* fix the elements over the diagonal */
		for (int i = 0; i < corr.length; i++) {
			for (int j = 1; j < corr[i].length; j++) {
				if (Math.abs(sigma[i] * sigma[i + j]) > Float.MIN_VALUE)
					corr[i][j] /= sigma[i] * sigma[i + j];
				
				assert(! Float.isNaN(corr[i][j]));
			}
		}

		sigma = null;

		/* fix the diagonal */
		for (int i = 0; i < corr.length; i++) {
			corr[i][0] = 1;
		}

		return corr;
	}
}
