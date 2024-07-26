package iemetrics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import GiciAnalysis.*;
import GiciException.*;
import GiciFile.*;

/**
 * @author Group on Interactive Coding of Images (GICI)
 * @version 1.0
 */
public class IEMetrics {
	/**
	 * Main method of IEMetrics application. It takes program arguments, loads images and compare them.
	 *
	 * @param args an array of strings that contains program parameters
	 */
	public static void main(String[] args)throws ErrorException{
		
//		Parse arguments
		ArgsParser parser = null;
		try{
			parser = new ArgsParser(args);
		}catch(Exception e){
			System.err.println("ARGUMENTS ERROR: " +  e.getMessage());
			e.printStackTrace();
			System.out.println("Please report this error (specifying image type and parameters) to: gici-dev@abra.uab.es");
			System.exit(1);
		}
		
		String inputFile1 = parser.getInputImageFile1();
		String inputFile2 = parser.getInputImageFile2();
		String resultsFolder = parser.getResultsFolder();
		
		int[] inputImageGeometry1 = parser.getInputImageGeometry1();
		int[] inputImageGeometry2 = parser.getInputImageGeometry2();
		
//		Image load
		
		LoadFile inputImage1 = null;
		LoadFile inputImage2 = null;
		
		try{
			if(LoadFile.isRaw(inputFile1)){				
				//Check parameters of image geometry
				if((inputImageGeometry1[0] <= 0) || (inputImageGeometry1[1] <= 0) || (inputImageGeometry1[2] <= 0)){
					throw new Exception("Image dimensions in \".raw\" or \".img\" data files must be positive (\"-h\" displays help).");
				}
				if((inputImageGeometry1[3] < 0) || (inputImageGeometry1[3] > 7)){
					throw new Exception("Image type in \".raw\" or \".img\" data must be between 0 to 7 (\"-h\" displays help).");
				}
				if((inputImageGeometry1[4] != 0) && (inputImageGeometry1[4] != 1)){
					throw new Exception("Image byte order  in \".raw\" or \".img\" data must be 0 or 1 (\"-h\" displays help).");
				}
				
				inputImage1 = new LoadFile(inputFile1, inputImageGeometry1[0], inputImageGeometry1[1], inputImageGeometry1[2], inputImageGeometry1[3], 
							inputImageGeometry1[4], false);				
			}else{
				inputImage1 = new LoadFile(inputFile1);	
			}
		}catch(Exception e){
			System.err.println("IMAGE LOADING ERROR: " + e.getMessage());
			System.exit(2);
		}
		
		try{
			if(LoadFile.isRaw(inputFile2)){				
				//Check parameters of image geometry
				if((inputImageGeometry2[0] <= 0) || (inputImageGeometry2[1] <= 0) || (inputImageGeometry2[2] <= 0)){
					throw new Exception("Image dimensions in \".raw\" or \".img\" data files must be positive (\"-h\" displays help).");
				}
				if((inputImageGeometry2[3] < 0) || (inputImageGeometry2[3] > 7)){
					throw new Exception("Image type in \".raw\" or \".img\" data must be between 0 to 7 (\"-h\" displays help).");
				}
				if((inputImageGeometry2[4] != 0) && (inputImageGeometry2[4] != 1)){
					throw new Exception("Image byte order  in \".raw\" or \".img\" data must be 0 or 1 (\"-h\" displays help).");
				}
				
				inputImage2 = new LoadFile(inputFile2, inputImageGeometry2[0], inputImageGeometry2[1], inputImageGeometry2[2], inputImageGeometry2[3], 
							inputImageGeometry2[4], false);				
			}else{
				inputImage2 = new LoadFile(inputFile2);	
			}
		}catch(Exception e){
			System.err.println("IMAGE LOADING ERROR: " + e.getMessage());
			System.exit(2);
		}
		
		float[][][] image1 = inputImage1.getImage();
		float[][][] image2 = inputImage2.getImage();
		
		// Check sizes
		int sizeZ = image1.length;
		int sizeY = image1[0].length;
		int sizeX = image1[0][0].length;
		
		if (sizeZ != image1.length || sizeY != image1[0].length || sizeX != image1[0][0].length) {
			System.err.println("IMAGE SIZES MUST BE EQUAL");
			System.exit(2);
		}
		
		try {
			// do something...
			
//			System.out.println("getMaximumSpectralSimilarity()sm.getMaximumSpectralAngle()"
//			+"sm.getMaximumSpectralInformationDivergence()sm.getMinimumPearsonsCorrelation()"
//			+"sm.getWangIndexLambda()sm.getWangIndexStack()sm.getWangIndexBoth()"
//			+"sm.getEskiciogluCubeFidelity()sm.getEskiciogluSpectralFidelity()"
//			+"sm.getEskiciogluStackFidelity()");
			
			if (true) {
				SpectralMetrics sm = (new SpectralMetrics(image1, image2)).call();

				System.out.print(
						sm.getMaximumSpectralSimilarity() + ":"
						+ sm.getMaximumSpectralAngle() + ":"
						+ sm.getMaximumSpectralInformationDivergence() + ":"
						+ sm.getMinimumPearsonsCorrelation() + ":"
						+ sm.getWangIndexLambda() + ":"
						+ sm.getWangIndexStack() + ":"
						+ sm.getWangIndexBoth() + ":"
						+ sm.getEskiciogluCubeFidelity() + ":"
						+ sm.getEskiciogluSpectralFidelity() + ":"
						+ sm.getEskiciogluStackFidelity()
				);

				sm = null;
			}
			/*
			 * POC preservation of classification
			 *  
			 * J. E. Fowler and J. T. Rucker, "3D Wavelet-Based Compression of Hyperspectral Imagery,"
			 * in Hyperspectral Data Exploitation: Theory and Applications, C.-I. Chang, Ed., chapter 14,
			 * pp. 379-407, John Wiley & Sons, Inc., Hoboken, NJ, 2007.
			 *
			 * Quote:
             * "We call the resulting distortion measure preservation of classification (POC),
             * which is measured as the percentage of pixels that do not change class due to compression."
             */

			/*
			 * A continuous extension could be: average(class_distance / max_distance)
			 * (TODO: think twice)
			 */
			
			// Lets go first with k-means/isodata
			
			float[][] points1 = new float[sizeX * sizeY][sizeZ];
			float[][] points2 = new float[sizeX * sizeY][sizeZ];
			
			for (int y = 0; y < sizeY; y++) {
				for (int x = 0; x < sizeX; x++) {
					for (int z = 0; z < sizeZ; z++) {
						points1[sizeX * y + x][z] = image1[z][y][x];
						points2[sizeX * y + x][z] = image2[z][y][x];
					}
				}
			}
			
			VectorDistance metricSAM = metricFactory.getSAM();
//			VectorDistance metricDot = metricFactory.getDot();
			VectorDistance metricEuler = metricFactory.getEuler();
			VectorDistance metricManhattan = metricFactory.getManhattan();
			
			float POC_K_MEANs_SAM = 0;
//			float POC_K_MEANs_Dot = 0;
			float POC_K_MEANs_Euler = 0;
			float POC_K_MEANs_Manhattan = 0;
			float POC_ISODATA_SAM = 0;
//			float POC_ISODATA_Dot = 0;
			float POC_ISODATA_Euler = 0;
			float POC_ISODATA_Manhattan = 0;
			
			// k-means options
			final int desiredClusters = 10;
			final int iterations = 100;
			
			// isodata options
			
			final int minimumSamplesPerCluster = 1;
			final int maximumPairToBeMerged = 2;
			
			final float deviationToSplit_SAM = 1;
			final float distanceToMerge_SAM = 0.5f;
//			final float deviationToSplit_Dot = 1;
//			final float distanceToMerge_Dot = 0.05f;
			final float deviationToSplit_Euler = 1;
			final float distanceToMerge_Euler = 5;
			final float deviationToSplit_Manhattan = 1;
			final float distanceToMerge_Manhattan = 5;
			
			testISODATA t;
			/* k-MEANs */
			if (true) {
				// 1
				t = new testISODATA (
						new IsoData(points1, metricSAM, desiredClusters, iterations),
						new IsoData(points2, metricSAM, desiredClusters, iterations));
				POC_K_MEANs_SAM = (t.call()).getPOC();

				if (resultsFolder != null)
					t.dumpFiles(resultsFolder, "poc_k_means_sam", sizeY, sizeX);

				// 2
//				t = new testISODATA (
//				new IsoData(points1, metricDot, desiredClusters, iterations),
//				new IsoData(points2, metricDot, desiredClusters, iterations));
//				POC_K_MEANs_Dot = (t.call()).getPOC();

//				if (resultsFolder != null)
//				t.dumpFiles(resultsFolder, "poc_k_means_dot", sizeY, sizeX);

				// 3
				t = new testISODATA (
						new IsoData(points1, metricEuler, desiredClusters, iterations),
						new IsoData(points2, metricEuler, desiredClusters, iterations));
				POC_K_MEANs_Euler = (t.call()).getPOC();

				if (resultsFolder != null)
					t.dumpFiles(resultsFolder, "poc_k_means_euler", sizeY, sizeX);

				// 4
				t = new testISODATA (
						new IsoData(points1, metricManhattan, desiredClusters, iterations),
						new IsoData(points2, metricManhattan, desiredClusters, iterations));
				POC_K_MEANs_Manhattan = (t.call()).getPOC();

				if (resultsFolder != null)
					t.dumpFiles(resultsFolder, "poc_k_means_manhattan", sizeY, sizeX);
			}
			
			/* ISODATAs */
			if (false) {
				// 1
				t = new testISODATA (
						new IsoData(points1, metricSAM, desiredClusters, iterations, minimumSamplesPerCluster, deviationToSplit_SAM, maximumPairToBeMerged, distanceToMerge_SAM),
						new IsoData(points2, metricSAM, desiredClusters, iterations, minimumSamplesPerCluster, deviationToSplit_SAM, maximumPairToBeMerged, distanceToMerge_SAM));
				POC_ISODATA_SAM = (t.call()).getPOC();

				if (resultsFolder != null)
					t.dumpFiles(resultsFolder, "poc_isodata_sam", sizeY, sizeX);

				// 2
//				t = new testISODATA (
//				new IsoData(points1, metricDot, desiredClusters, iterations, minimumSamplesPerCluster, deviationToSplit_Dot, maximumPairToBeMerged, distanceToMerge_Dot),
//				new IsoData(points2, metricDot, desiredClusters, iterations, minimumSamplesPerCluster, deviationToSplit_Dot, maximumPairToBeMerged, distanceToMerge_Dot));
//				POC_ISODATA_Dot = (t.call()).getPOC();

//				if (resultsFolder != null)
//				t.dumpFiles(resultsFolder, "poc_isodata_dot", sizeY, sizeX);

				// 3
				t = new testISODATA (
						new IsoData(points1, metricEuler, desiredClusters, iterations, minimumSamplesPerCluster, deviationToSplit_Euler, maximumPairToBeMerged, distanceToMerge_Euler),
						new IsoData(points2, metricEuler, desiredClusters, iterations, minimumSamplesPerCluster, deviationToSplit_Euler, maximumPairToBeMerged, distanceToMerge_Euler));
				POC_ISODATA_Euler = (t.call()).getPOC();

				if (resultsFolder != null)
					t.dumpFiles(resultsFolder, "poc_isodata_euler", sizeY, sizeX);

				// 4
				t = new testISODATA (
						new IsoData(points1, metricManhattan, desiredClusters, iterations, minimumSamplesPerCluster, deviationToSplit_Manhattan, maximumPairToBeMerged, distanceToMerge_Manhattan),
						new IsoData(points2, metricManhattan, desiredClusters, iterations, minimumSamplesPerCluster, deviationToSplit_Manhattan, maximumPairToBeMerged, distanceToMerge_Manhattan));
				POC_ISODATA_Manhattan = (t.call()).getPOC();

				if (resultsFolder != null)
					t.dumpFiles(resultsFolder, "poc_isodata_manhattan", sizeY, sizeX);
			}
			
			// Free some memory			
			t = null;
			points1 = null;
			points2 = null;
			
			System.out.print(":" 
					+ POC_K_MEANs_SAM + ":"
					//+ POC_K_MEANs_Dot + ":"
					+ POC_K_MEANs_Euler + ":"
					+ POC_K_MEANs_Manhattan + ":"
					+ POC_ISODATA_SAM + ":"
					//+ POC_ISODATA_Dot + ":"
					+ POC_ISODATA_Euler + ":"
					+ POC_ISODATA_Manhattan);
			
			// RX classification
			
			float[][] rx1 = (new RXAlgorithm(image1)).call();
			float[][] rx2 = (new RXAlgorithm(image2)).call();

			if (resultsFolder != null) {
				float[][][] rx1p = {rx1};
				float[][][] rx2p = {rx2};
				
				SaveFile.SaveFileFormat(rx1p, resultsFolder + "/rx1.png", 2);
				SaveFile.SaveFileFormat(rx2p, resultsFolder + "/rx2.png", 2);
			}
			
			// Find threshold
			/* Automated method from:
			 * Anomaly Detection and Classification for Hyperspectral Imagery
			 * Chein-I Chang, Senior Member, IEEE, and Shao-Shan Chiang, Member, IEEE
			 */
			
			final List<Float> rxThresholdList1 = new ArrayList<Float>();
			final List<Float> rxThresholdList2 = new ArrayList<Float>();
			
			for (int y = 0; y < rx1.length; y++) {
				for (int x = 0; x < rx1[0].length; x++) {
					rxThresholdList1.add(rx1[y][x]);
					rxThresholdList2.add(rx2[y][x]);
				}
			}
			
			Collections.sort(rxThresholdList1);
			Collections.sort(rxThresholdList2);
			
			final float confidenceCoefficient = 0.99f;
			final float rxThreshold1 = rxThresholdList1.get((int)Math.floor((rxThresholdList1.size() - 1) * confidenceCoefficient));
			final float rxThreshold2 = rxThresholdList2.get((int)Math.floor((rxThresholdList2.size() - 1) * confidenceCoefficient));
			
			int rxChanged = 0;
			int rxTotal = rx1.length * rx1[0].length;
			
			for (int y = 0; y < rx1.length; y++) {
				for (int x = 0; x < rx1[0].length; x++) {
					if ((rx1[y][x] < rxThreshold1 && rx2[y][x] >= rxThreshold2)
							|| (rx1[y][x] >= rxThreshold1 && rx2[y][x] < rxThreshold2)
							|| (Float.isNaN(rxThreshold1) && ! Float.isNaN(rxThreshold2))
							|| (! Float.isNaN(rxThreshold1) && Float.isNaN(rxThreshold2)) ) {
						rxChanged++;
					}
				}
			}
			
			System.out.print(":" + ((rxTotal - rxChanged) / (double)rxTotal));
			
			System.out.println("");
		}catch(Exception e){
			e.printStackTrace();
			System.err.println("PROCESSING ERROR: " + e.getMessage());
			System.exit(2);
		}
	}
}
