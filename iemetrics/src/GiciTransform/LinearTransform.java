package GiciTransform;

import GiciMatrix.MatrixAlgebra;

/**
 * @author Ian Blanes
 *
 */
public class LinearTransform {
	/**
	 * Transform that this object will apply.
	 */
	private float[][] transform;

	/**
	 * The dimension that this object will work on.
	 */
	private int dimension;

	/**
	 * Constructor.
	 *
	 * @param t Transform Matrix.
	 * @param d Dimension.
	 */
	public LinearTransform(final float[][] t, final int d) {
		this.transform = t;
		this.dimension = d;

		final int maxDimensions = 3;

		assert (t.length > 0 && t.length == t[0].length);
		assert (dimension < maxDimensions && dimension >= 0);
	}

	/**
	 * Applies the transform in the same input.
	 *
	 * @param imageSamples input and output where the transform is performed.
	 */
	public final void transformInPlace(final float[][][] imageSamples) {

		int d0 = imageSamples.length;
		int d1 = imageSamples[0].length;
		int d2 = imageSamples[0][0].length;
		float[] v;

		switch (dimension) {
		case 0:
			assert (d0 == transform.length);
			v = new float[d0];

			for (int i = 0; i < d1; i++) {
				for (int j = 0; j < d2; j++) {
					for (int k = 0; k < d0; k++) {
						v[k] = imageSamples[k][i][j];
					}

					v = MatrixAlgebra.multiplicationCV(transform, v);

					for (int k = 0; k < d0; k++) {
						imageSamples[k][i][j] = v[k];
					}
				}
			}
			break;
		case 1:
			assert (d1 == transform.length);
			v = new float[d1];

			for (int i = 0; i < d0; i++) {
				for (int j = 0; j < d2; j++) {
					for (int k = 0; k < d1; k++) {
						v[k] = imageSamples[i][k][j];
					}

					v = MatrixAlgebra.multiplicationCV(transform, v);

					for (int k = 0; k < d1; k++) {
						imageSamples[i][k][j] = v[k];
					}
				}
			}
			break;
		case 2:
			assert (d2 == transform.length);
			v = new float[d2];

			for (int i = 0; i < d0; i++) {
				for (int j = 0; j < d1; j++) {
					for (int k = 0; k < d2; k++) {
						v[k] = imageSamples[i][j][k];
					}

					v = MatrixAlgebra.multiplicationCV(transform, v);

					for (int k = 0; k < d2; k++) {
						imageSamples[i][j][k] = v[k];
					}
				}
			}
			break;
		default:
			throw new ArrayIndexOutOfBoundsException();
		}
	}

	/**
	 * Applies the transform and returns the result.
	 *
	 * @param imageSamples input where the transform will be computed.
	 * @return result.
	 */
	public final float[][][] transform(final float[][][] imageSamples) {
		float[][][] r = new float[imageSamples.length][][];

		for (int z = 0; z < imageSamples.length; z++) {
			r[z] = new float[imageSamples[z].length][];

			for (int x = 0; x < imageSamples[z].length; x++) {
				r[z][x] = new float[imageSamples[z][x].length];

				System.arraycopy(imageSamples[z][x], 0, r[z][x], 0, imageSamples[z][x].length);
			}
		}

		transformInPlace(r);

		return r;
	}
}
