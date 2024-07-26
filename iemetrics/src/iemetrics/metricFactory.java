package iemetrics;

import GiciAnalysis.VectorDistance;

public class metricFactory {
	static VectorDistance getSAM () {
		return new VectorDistance(){
			public float distance(float[] a, float[] b) {
				/* dot product */
				float acc1 = 0;
				float acc2 = 0;
				float acc3 = 0;

				for (int i = 0; i < a.length; i++) {
					acc1 += a[i] * b[i];
					acc2 += a[i] * a[i];
					acc3 += b[i] * b[i];
				}

				// FIXME!! a metric has to be defined for every point in space, and
				// this one is not!
				if (acc2 * acc3 < 0.000001f) {
					return 0;
				}
				
				double t = acc1 / (Math.sqrt(acc2) * Math.sqrt(acc3)); 
				
				// Check for slightly larger dot products that produce NaNs on acos.
				if (Math.abs(t) > 1) {
					return (float) Math.acos(Math.signum(t));
				}
				
				return (float) Math.acos(t);
				// This result is in degrees!
			}
		};
	}

	static VectorDistance getDot() {
		return new VectorDistance(){
			public float distance(float[] a, float[] b) {
				/* dot product */
				float acc = 0;

				for (int i = 0; i < a.length; i++) {
					acc += a[i] * b[i];
				}

				return acc;
			}
		};
	}

	static VectorDistance getEuler() {
		return new VectorDistance(){
			public float distance(float[] a, float[] b) {
				float acc = 0;

				for (int i = 0; i < a.length; i++) {
					float d = a[i] - b[i];
					acc += d * d;
				}

				return (float)Math.sqrt(acc);
			}
		};
	}

	static VectorDistance getManhattan() {
		return new VectorDistance(){
			public float distance(float[] a, float[] b) {
				float acc = 0;

				for (int i = 0; i < a.length; i++) {
					acc += Math.abs(a[i] - b[i]);
				}

				return acc;
			}
		};
	}
}
