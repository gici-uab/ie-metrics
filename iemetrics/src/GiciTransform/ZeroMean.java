package GiciTransform;

public class ZeroMean {
	
	private int dimension;
	private boolean forward;
	private boolean meansCorrect = false;
	private boolean reversible;
	
	private final float[][][] image;
	private float[] means;
	
	/**
	 * 
	 * @param image
	 * @param dimension 0z .. 2x
	 * @param reversible
	 */
	
	public ZeroMean(final float[][][] image, int dimension, boolean reversible) {
		this.image = image;
		this.dimension = dimension;
		this.forward = true;
		this.reversible = reversible;
		
		assert(dimension >= 0 && dimension <= 2);
	}
	
	public ZeroMean(float[][][] image, int dimension, float[] means) {
		this.image = image;
		this.dimension = dimension;
		this.forward = false;
		this.means = means;
		this.reversible = true;
		
		assert(dimension >= 0 && dimension <= 2);
	}
	
	private void calculateMeans () {
		
		if (meansCorrect) { return; }
		
		meansCorrect = true;
		
		// ! FWD
		
		if (! forward) {
			// Reverse
			for (int i = 0; i < means.length; i++) {
				means[i] = -means[i]; 
			}

			return;
		}
		
		// FWD
		
		double[] dmeans = null;

		switch (dimension) {
		case 0:
			dmeans = new double[image.length];
			means = new float[image.length];

			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						dmeans[i] += image[i][j][k];
					}
				}
			}

			for (int i = 0; i < dmeans.length; i++) {
				means[i] = (float) (dmeans[i] / (double) (image[0].length * image[0][0].length));
			}				

			break;
		case 1:
			dmeans = new double[image[0].length];
			means = new float[image[0].length];

			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						dmeans[j] += image[i][j][k];
					}
				}
			}

			for (int i = 0; i < dmeans.length; i++) {
				means[i] = (float) (dmeans[i] / (double) (image.length * image[0][0].length));
			}

			break;
		case 2:
			dmeans = new double[image[0][0].length];
			means = new float[image[0][0].length];

			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						dmeans[k] += image[i][j][k];
					}
				}
			}

			for (int i = 0; i < dmeans.length; i++) {
				means[i] = (float) (dmeans[i] / (double) (image.length * image[0].length));
			}

			break;
		}

		if (reversible) {			
			for (int i = 0; i < dmeans.length; i++) {
				means[i] = Math.round(means[i]);
			}	
		}
	}
	
	public float[][][] applyMeans () {
		calculateMeans();
		
		// Apply
		float[][][] result = new float[image.length][image[0].length][image[0][0].length];
		
		switch (dimension) {
		case 0:
			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						result[i][j][k] = image[i][j][k] - means[i];
					}
				}
			}
			break;
		case 1:
			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						result[i][j][k] = image[i][j][k] - means[j];
					}
				}
			}
			break;
		case 2:
			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						result[i][j][k] = image[i][j][k] - means[k];
					}
				}
			}
			break;
		}
		
		return result;
	}
	
	public void applyMeansInPlace () {
		calculateMeans();
		
		// Apply
		
		switch (dimension) {
		case 0:
			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						image[i][j][k] -= means[i];
					}
				}
			}
			break;
		case 1:
			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						image[i][j][k] -= means[j];
					}
				}
			}
			break;
		case 2:
			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[i].length; j++) {
					for (int k = 0; k < image[i][j].length; k++) {
						image[i][j][k] -= means[k];
					}
				}
			}
			break;
		}
		
	}

	public float[] getMeans() {
		assert(forward);

		calculateMeans();
		
		return means;
	}
}
