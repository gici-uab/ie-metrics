package GiciAnalysis;

import java.util.concurrent.Callable;

import GiciMatrix.MatrixAlgebra;
import GiciTransform.KarhunenLoeveTransform;
import GiciTransform.LinearTransform;
import GiciTransform.ZeroMean;

// Not optimized at all, so it might be highly expensive

public class RXAlgorithm implements Callable {
	// Source Data
	final float[][][] image;
	final int sizeX;
	final int sizeY;
	final int sizeZ;
	
	public RXAlgorithm(final float[][][] image) {
		super();

		this.image = image;

		this.sizeZ = image.length;
		this.sizeY = image[0].length;
		this.sizeX = image[0][0].length;
	}

	public float[][] call() throws Exception {
		// Means (FIXME: either optimize or generalize)
		float[][][] means = new float[sizeZ][sizeY][sizeX]; 
		
		if (false) {
			// a) local window
			int w = 3;
			
			for (int y = 0; y < sizeY; y++) {
				for (int x = 0; x < sizeX; x++) {
					float[] v = new float[sizeZ];

					int sy = Math.max(0, y - w);
					int sx = Math.max(0, x - w);
					int ey = Math.min(sizeY, y + w + 1);
					int ex = Math.min(sizeX, x + w + 1);

					for (int dy = sy; dy < ey; dy++) {
						for (int dx = sx; dx < ex; dx++) {
							for (int z = 0; z < sizeZ; z++) {
								v[z] += image[z][dy][dx];
							}						
						}
					}

					for (int z = 0; z < sizeZ; z++) {
						means[z][y][x] = v[z] / ((ey - sy) * (ex - sx)); 
					}
				}
			}
		} else {
			// b) whole image
			
			float[] v = new float[sizeZ];
			
			for (int z = 0; z < sizeZ; z++) {
				for (int y = 0; y < sizeY; y++) {
					for (int x = 0; x < sizeX; x++) {
						v[z] += image[z][y][x];
					}
				}
			}
			
			for (int z = 0; z < sizeZ; z++) {
				v[z] /= sizeX * sizeY; 
			}	
			
			for (int z = 0; z < sizeZ; z++) {
				for (int y = 0; y < sizeY; y++) {
					for (int x = 0; x < sizeX; x++) {
						means[z][y][x] = v[z];
					}
				}
			}
		}
		
		// Substract the global image mean
		float[][][] imagePerPointZeroMean = new float[sizeZ][sizeY][];

		for (int z = 0; z < sizeZ; z++) {
			for (int y = 0; y < sizeY; y++) {
				// save some memory
				imagePerPointZeroMean[z][y] = new float[sizeX];
				
				for (int x = 0; x < sizeX; x++) {
					imagePerPointZeroMean[z][y][x] = image[z][y][x] - means[z][y][x];
				}
				
				// save some memory
				means[z][y] = null;
			}
		}			

		means = null;
		
		// Covariance Matrix
		final float[][] cov = ImageCovariance.generateCovarianceMatrix(image, 0);
		
		// klt space
		float[][][] klt = KarhunenLoeveTransform.generateTransform(cov);
		LinearTransform lt = new LinearTransform(MatrixAlgebra.transposeC(klt[0]),0);
		lt.transformInPlace(imagePerPointZeroMean);

		// Eulerian distance
		float[][] distance = new float[sizeY][sizeX];

		for (int i = 0; i < sizeY; i++) {
			for (int j = 0; j < sizeX; j++) {
				float acc = 0;

				for (int k = 0; k < sizeZ; k++) {
					float d = (imagePerPointZeroMean[k][i][j]);
					float r = d * d / klt[1][k][k];
					
					// Only if d * d = 0 && klt = 0 skip this result
					// d * d > 0 && klt = 0 might still be valid (?)
					if (Float.isNaN(r))
						continue;
					
					acc += r;
				}

				distance[i][j] = (float) Math.sqrt(acc);
				
				assert(! Float.isNaN(distance[i][j]));
			}
		}
		
		return distance;
	}

	public float[][] callold() throws Exception {
		// Per col zero mean
		ZeroMean zm = new ZeroMean(image, 0, false);
		final float[][][] imagePerColZeroMean = zm.applyMeans();
		
		// Result
		float[][] distance;
		
		// Means (FIXME: either optimize or generalize)
		float[][][] means = new float[sizeZ][sizeY][sizeX]; 
		
		if (false) {
			// a) local window
			int w = 3;
			
			for (int y = 0; y < sizeY; y++) {
				for (int x = 0; x < sizeX; x++) {
					float[] v = new float[sizeZ];

					int sy = Math.max(0, y - w);
					int sx = Math.max(0, x - w);
					int ey = Math.min(sizeY, y + w + 1);
					int ex = Math.min(sizeX, x + w + 1);

					for (int dy = sy; dy < ey; dy++) {
						for (int dx = sx; dx < ex; dx++) {
							for (int z = 0; z < sizeZ; z++) {
								v[z] += image[z][dy][dx];
							}						
						}
					}

					for (int z = 0; z < sizeZ; z++) {
						means[z][y][x] = v[z] / ((ey - sy) * (ex - sx)); 
					}
				}
			}
		} else {
			// b) whole image
			
			float[] v = new float[sizeZ];
			
			for (int z = 0; z < sizeZ; z++) {
				for (int y = 0; y < sizeY; y++) {
					for (int x = 0; x < sizeX; x++) {
						v[z] += image[z][y][x];
					}
				}
			}
			
			for (int z = 0; z < sizeZ; z++) {
				v[z] /= sizeX * sizeY; 
			}	
			
			for (int z = 0; z < sizeZ; z++) {
				for (int y = 0; y < sizeY; y++) {
					for (int x = 0; x < sizeX; x++) {
						means[z][y][x] = v[z];
					}
				}
			}
		}
		
		// Covariance
		float[][] cov = KarhunenLoeveTransform.toCompleteMirrorFill(ImageCovariance.generateCovarianceMatrix(imagePerColZeroMean, 0));
		//MatrixAlgebra.printMatrix(cov);
		
		assert(Math.abs(MatrixAlgebra.determinant(cov)) > 0.0001f);
		
		float[][] cov1 = MatrixAlgebra.inverseC(cov);
		//MatrixAlgebra.printMatrix(cov1);
		distance = new float[sizeY][sizeX];
		
		for (int y = 0; y < sizeY; y++) {
			for (int x = 0; x < sizeX; x++) {
				float[][] v = new float[sizeZ][1];
				
				for (int z = 0; z < sizeZ; z++) {
					v[z][0] = image[z][y][x] - means[z][y][x]; 
				}
				
				float[][] vt = MatrixAlgebra.transposeC(v);
				float[][] r = MatrixAlgebra.multiplicationCC(vt, MatrixAlgebra.multiplicationCC(cov1, v));
				
				assert(r.length == 1 && r[0].length == 1);
				
				distance[y][x] = (float)Math.sqrt(r[0][0]);
			}
		}
		
		return distance;
	}


	
}
