package iemetrics;

/**
 * Arguments parser for ffc. This class analyses a string of arguments and extract and check its validity.
 * Usage example:<br>
 * &nbsp; construct<br>
 * &nbsp; [showArgsInfo]<br>
 * &nbsp; [get functions]<br>
 *
 * @author Group on Interactive Coding of Images (GICI)
 * @version 1.0
 */
public class ArgsParser{

	/**
	 * Arguments specificiation. The array describes argument, explain what is used and its default parameters. First index of array is argument; second specifies:<br>
	 *   <ul>
	 *     <li> 0 - short argument specification (i.e. "-i")
	 *     <li> 1 - long argument specification (i.e. "--inputFile")
	 *     <li> 2 - parsing specification of argument ({} indicates mandatority, [] optionality)
	 *     <li> 3 - default values
	 *     <li> 4 - mandatory argument ("1") or non mandatory argument ("0")
	 *     <li> 5 - explanation
	 *   </ul>
	 * <p>
	 * String arguments.
	 */
	String[][] argsSpecification = {
		{"-h", "--help", "", "", "0",
			"Displays this help and exits program."
		},
		{"-i1", "--inputImageFile1", "{string}", "", "1",
			"Input image. Valid formats are: pgm, ppm, pbm, jpg, tiff, png, bmp, gif, fpx. If image is raw data file extension must be \".raw\" or \".img\" and \"-ig\" parameter is mandatory."
		},
		{"-ig1", "--inputImageGeometry1", "{int int int int boolean}", "", "0",
			"Geometry of input raw image data. Parameters are:\n    1- zSize (number of image components)\n    2- ySize (image height)\n    3- xSize (image width)\n    4- data type. Possible values are:\n \t 0- boolean (1 byte)\n \t 1- unsigned int (1 byte)\n \t 2- unsigned int (2 bytes)\n \t 3- signed int (2 bytes)\n \t 4- signed int (4 bytes)\n \t 5- signed int (8 bytes)\n \t 6- float (4 bytes)\n \t 7- double (8 bytes)\n    5- Byte order (0 if BIG ENDIAN, 1 if LITTLE ENDIAN)."
		},
		{"-i2", "--inputImageFile2", "{string}", "", "1",
			"Output image file name. Valid formats are: pgm, ppm, pbm, jpg, tiff, png, bmp, gif, fpx. If image is raw data file extension must be \".raw\" or \".img\" and \"-og\" parameter is mandatory."
		},
		{"-ig2", "--inputImageGeometry2", "{int int int int boolean}", "", "0",
			"Geometry of output raw image data. Parameters are:\n    1- zSize (number of image components)\n    2- ySize (image height)\n    3- xSize (image width)\n    4- data type. Possible values are:\n \t 0- boolean (1 byte)\n \t 1- unsigned int (1 byte)\n \t 2- unsigned int (2 bytes)\n \t 3- signed int (2 bytes)\n \t 4- signed int (4 bytes)\n \t 5- signed int (8 bytes)\n \t 6- float (4 bytes)\n \t 7- double (8 bytes)\n    5- Byte order (0 if BIG ENDIAN, 1 if LITTLE ENDIAN)."
		},
		{"-dr", "--dumpResults", "{string}", "", "0",
			"Folder on which to dump various results."
		},
	};

	//ARGUMENTS VARIABLES
	String inputImageFile1 = null;
	String inputImageFile2 = null;
	int[] inputImageGeometry1 = null;
	int[] inputImageGeometry2 = null;
	String resultsFolder = null;
	
	 /**
	  * Class constructor that receives the arguments string and initializes all the arguments
	  * 
	  * @param args the array of strings passed at the command line
	  * 
	  * @throws Exception when an invalid parsing is detected or some problem with method invocation occurs
	  */
	public ArgsParser(String[] args) throws Exception{
		int argNum = 0;
		boolean[] argsFound = new boolean[argsSpecification.length];

		//Arguments parsing
		for(int i = 0; i < argsSpecification.length; i++){
			argsFound[i] = false;
		}
		while(argNum < args.length){
			int argFound = argFind(args[argNum]);
			if(argFound != -1){
				if(!argsFound[argFound]){
					argsFound[argFound] = true;
					int argOptions = argNum + 1;
					while(argOptions < args.length){
						if(argFind(args[argOptions]) != -1){
							break;
						}else{
							argOptions++;
						}
					}
					int numOptions = argOptions - argNum;
					String[] options = new String[numOptions];
					System.arraycopy(args, argNum, options, 0, numOptions);
					argNum = argOptions;
					switch(argFound){
					case 0: //-h  --help
						showArgsInfo();
						System.exit(1);
						break;
					case 1: //-i1  --inputImageFile1
						inputImageFile1 = parseString(options);
						if(inputImageFile1.endsWith(".raw")){
							argsSpecification[2][4] = "1";
						}
						break;
					case 2: //-ig1  --inputImageGeometry1
						inputImageGeometry1 = parseIntegerArray(options, 5);
						break;
					case 3: //-i2  --inputImageFile2
						inputImageFile2 = parseString(options);
						if(inputImageFile2.endsWith(".raw")){
							argsSpecification[4][4] = "1";
						}
						break;
					case 4: //-ig2  --inputImageGeometry2
						inputImageGeometry2 = parseIntegerArray(options, 5);
						break;
					case 5: //-dr --dumpResults
						resultsFolder = parseString(options);
						break;	
					}
				}else{
					throw new Exception("Argument \"" + args[argNum] + "\" repeated.");
				}
			}else{
				throw new Exception("Argument \"" + args[argNum] + "\" unrecognized.");
			}
		}

		//Check mandatory arguments
		for(int i = 0; i < argsSpecification.length; i++){
			if(argsSpecification[i][4].compareTo("1") == 0){
				if(!argsFound[i]){
					throw new Exception("Argument \"" + argsSpecification[i][0] + "\" is mandatory (\"-h\" displays help).");
				}
			}
		}
	}

	/**
	 * Finds the argument string in arguments specification array.
	 *
	 * @param arg argument to find out in argsSpecification
	 * @return the argument index of argsSpecification (-1 if it doesn't exist)
	 */
	int argFind(String arg){
		int argFound = 0;
		boolean found = false;

		while((argFound < argsSpecification.length) && !found){
			if((arg.compareTo(argsSpecification[argFound][0]) == 0) || (arg.compareTo(argsSpecification[argFound][1]) == 0)){
				found = true;
			}else{
				argFound++;
			}
		}
		return(found ? argFound: -1);
	}

	/**
	 * This function shows arguments information to console.
	 */
	public void showArgsInfo(){
		System.out.println("Arguments specification: ");
		for(int numArg = 0; numArg < argsSpecification.length; numArg++){
			char beginMandatory = '{', endMandatory = '}';
			if(argsSpecification[numArg][4].compareTo("0") == 0){
				//No mandatory argument
				beginMandatory = '[';
				endMandatory = ']';
			}
			System.out.print("\n" + beginMandatory + " ");
			System.out.print("{" + argsSpecification[numArg][0] + "|" + argsSpecification[numArg][1] + "} " + argsSpecification[numArg][2]);
			System.out.println(" " + endMandatory);
			System.out.println("  Explanation:\n    " + argsSpecification[numArg][5]);
			System.out.println("  Default value: " + argsSpecification[numArg][3]);
		}
	}


	/////////////////////
	//PARSING FUNCTIONS//
	/////////////////////
	//These functions receives a string array that contains in first position the argument and then their options//

	int parseIntegerPositive(String[] options) throws Exception{
		int value = 0;

		if(options.length == 2){
			try{
				value = Integer.parseInt(options[1]);
				if(value < 0){
					throw new Exception("\"" + options[1] + "\" of argument \"" + options[0] + "\" is must be a positive integer.");
				}
			}catch(NumberFormatException e){
				throw new Exception("\"" + options[1] + "\" of argument \"" + options[0] + "\" is not a parsable integer.");
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes one option. Try \"-h\" to display help.");
		}
		return(value);
	}

	String parseString(String[] options) throws Exception{
		String value = "";

		if(options.length == 2){
			value = options[1];
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes one option. Try \"-h\" to display help.");
		}
		return(value);
	}

	int[] parseIntegerArray(String[] options) throws Exception{
		int[] value = null;

		if(options.length >= 2){
			value = new int[options.length - 1];
			for(int numOption = 1; numOption < options.length; numOption++){
				try{
						value[numOption - 1] = Integer.parseInt(options[numOption]);
				}catch(NumberFormatException e){
					throw new Exception("\"" + options[numOption] + "\" of argument \"" + options[0] + "\" is not a parsable integer.");
				}
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes one or more options. Try \"-h\" to display help.");
		}
		return(value);
	}

	int[] parseIntegerArray(String[] options, int numOptions) throws Exception{
		int[] value = null;

		if(options.length == numOptions+1){
			value = new int[options.length - 1];
			for(int numOption = 1; numOption < options.length; numOption++){
				try{
						value[numOption - 1] = Integer.parseInt(options[numOption]);
				}catch(NumberFormatException e){
					throw new Exception("\"" + options[numOption] + "\" of argument \"" + options[0] + "\" is not a parsable integer.");
				}
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes " + numOptions +" options. Try \"-h\" to display help.");
		}
		return(value);
	}
	
	float[] parseFloatArray(String[] options) throws Exception{
		float[] value = null;

		if(options.length >= 2){
			value = new float[options.length - 1];
			for(int numOption = 1; numOption < options.length; numOption++){
				try{
						value[numOption - 1] = Float.parseFloat(options[numOption]);
				}catch(NumberFormatException e){
					throw new Exception("\"" + options[numOption] + "\" of argument \"" + options[0] + "\" is not a parsable float.");
				}
			}
		}else{
			throw new Exception("Argument \"" + options[0] + "\" takes one or more options. Try \"-h\" to display help.");
		}
		return(value);
	}


	///////////////////////////
	//ARGUMENTS GET FUNCTIONS//
	///////////////////////////

	public String getInputImageFile1(){
		return(inputImageFile1);
	}
	public String getInputImageFile2(){
		return(inputImageFile2);
	}
	public int[] getInputImageGeometry1(){
		return(inputImageGeometry1);
	}
	public int[] getInputImageGeometry2(){
		return(inputImageGeometry2);
	}
	public String getResultsFolder(){
		return(resultsFolder);
	}
}

