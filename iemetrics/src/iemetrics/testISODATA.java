package iemetrics;

import java.util.concurrent.Callable;

import GiciAnalysis.ClassToColor;
import GiciAnalysis.IsoData;
import GiciAnalysis.PreservationOfClassification;
import GiciFile.SaveFile;

public class testISODATA implements Callable {
	boolean called = false;
	
	final IsoData test1, test2;

	int[] class1;
	int[] class2;
	float POC;
	
	public testISODATA (IsoData _test1, IsoData _test2) {
		test1 = _test1;
		test2 = _test2;
	}

	public testISODATA call() throws Exception {
		called = true;
		
		class1 = test1.call();
		class2 = test2.call();
		
		POC = (new PreservationOfClassification(class1, class2)).call();
		
		return this;
	}
	
	public float getPOC() {
		assert(called);
		
		return POC;
	}
	
	public void dumpFiles(String folder, String prefix, int y, int x) throws Exception {
		assert (called);
		
		ClassToColor cc = new ClassToColor(class1, y, x, 255);
		float[][][] i1 = cc.call();
		
		cc = new ClassToColor(class2, y, x, 255);
		float[][][] i2 = cc.call();
		
		SaveFile.SaveFileFormat(i1, folder + "/" + prefix + "1.png", 2);
		SaveFile.SaveFileFormat(i2, folder + "/" + prefix + "2.png", 2);
	}
}
