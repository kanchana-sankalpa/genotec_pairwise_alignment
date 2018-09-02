/**
 *  Local Alignment with Modified Algorithm
 *  @Genotec
 */


package alignment_test;


import java.io.BufferedReader;
import java.io.FileReader;

//import com.genotec.util.Constants;

public class GT_LocalAlignmentManager {

	private String refSeq = "";
	private String patientSeq1 = "";
	private String patientSeq2 = "";

	private int[][] scoreMatrix = null;
	private String[] pathMatrix = null;
	public static String ref1 = null;
	public static String target1 = null;
	public static String target2 = null;
	public static int startIndex = 0;
	public static int endIndex = 0;

	final int GAP_PANELTY = -1;
	
	public GT_LocalAlignmentManager (String patientSeq1, String patientSeq2,String diseaseId){
		this.patientSeq1 = patientSeq1;
		this.patientSeq2 = patientSeq2;
		
		//start test data
		/*this.patientSeq1 = "TCATACATTTTTCTCTAACTGCAAACATAATGTTTTCCCTTGTATTTTACAGATGCAAAG"
			+ "AGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTAA"
			+ "CATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATC"
			+ "CTTCCTTGGTAAAACCATTTGTTTTCTTCTTCTTCTTCTTCTTCTTTTCTTTTTTTTTTCTT"
			+ "TTTTTTTTTTGAGATGGAGTCTTGCTCTGTGGCCCAGGCTAGAAGCAGTCCTCCTGCCTTAG"
			+ "CCCCCTTAGTAGCTGGGATTACAGGCACGCGCCACCATGCCAGGCTAATTTTTGTATTTTTA"
			+ "GTAGAGACGGGGTTTCATCATGTTGGCCAGGCTGGTCTCGAACTCCTAACCTCAGGTGATCC"
			+ "ACCCACCTCGGCTCCCCAAATTGCTGGGATTACAGGTGTGAGCCACTGTGCCCGGCCGGTAA"
			+ "AACCATTTTCATTTATTCTGGCAACATCTCTTTATTGAGCATTGTGAATATGTTAGTGAATG"
			+ "TGCTAGATGCTCATAGATTTATATAAAAAGTTAGTGAAGAAGGAAAGATGGTATATTAAGTG"
			+ "GTTAGACAAGTGTTCTAATCAGTTAGAGTTCAGAGAAGGTCAGGGTACCTGATATAATCAAG"
			+ "AGAGAGACCTTACAGCCAGGTGAGGTGAATGTACCTATAATCCCAGCTACTTAGGAGGCT";
		this.patientSeq2 = "TCATACATTTTTCTCTAACTGCAAACATAATGTTTTCCCTTGTATTTTACAGATGCAAAC"
		+ "AGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTAT"
		+ "CATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATC"
		+ "CTTCCTTGGTAAAACCATTTGTTTTCTTCTTCTTCTTCTTCTTCTTTTCTTTTTTTTTTCTT"
		+ "TTTTTTTTTTGAGATGGAGTCTTGCTCTGTGGCCCAGGCTAGAAGCAGTCCTCCTGCCTTAG"
		+ "CCCCCTTAGTAGCTGGGATTACAGGCACGCGCCACCATGCCAGGCTAATTTTTGTATTTTTA"
		+ "GTAGAGACGGGGTTTCATCATGTTGGCCAGGCTGGTCTCGAACTCCTAACCTCAGGTGATCC"
		+ "ACCCACCTCGGCTCCCCAAATTGCTGGGATTACAGGTGTGAGCCACTGTGCCCGGCCGGTAA"
		+ "AACCATTTTCATTTATTCTGGCAACATCTCTTTATTGAGCATTGTGAATATGTTAGTGAATG"
		+ "TGCTAGATGCTCATAGATTTATATAAAAAGTTAGTGAAGAAGGAAAGATGGTATATTAAGTG"
		+ "GTTAGACAAGTGTTCTAATCAGTTAGAGTTCAGAGAAGGTCAGGGTACCTGATATAATCAAG"
		+ "AGAGAGACCTTACAGCCAGGTGAGGTGAATGTACCTATAATCCCAGCTACTTAGGAGGCT";*/
		//end test data
		
		ref1 = null;
		target1 = null;
		target2 = null;
		startIndex = 0;
		this.localAlignment(diseaseId);
	}
	
	/**
	 * 
	 * @param diseaseId
	 */
	private void localAlignment(String diseaseId) {
		this.refSeq = readReferenceFile(diseaseId);
		char[] refArr = this.refSeq.toCharArray();
		char[] patientArr1 = this.patientSeq1.toCharArray();
		char[] patientArr2 = this.patientSeq2.toCharArray();
		scoreMatrix = new int[2][patientArr1.length + 1];//scoring matrix (2*patient sequence length+1)
		pathMatrix = new String[refArr.length];// store path to use in trace back
		initialize(refArr, patientArr1, patientArr2);//initialize the scoring
	}

	/**
	 * Read the disease reference file from repository
	 * @param diseaseId
	 * @return
	 */
	private String readReferenceFile(String diseaseId) {
		String querySeq = null;
		try {
			BufferedReader br = new BufferedReader(new FileReader(diseaseId));
			String line = br.readLine();
			while (line != null) {
				if (line != null) {
					querySeq += line;
				}
				line = br.readLine();
			}
			br.close();
			System.out.println("File readed...");

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return querySeq;
	}



	/**
	 * Initialize the scoring process
	 * @param refArr
	 * @param patientArr1
	 * @param patientArr2
	 */
	public void initialize(char[] refArr, char[] patientArr1, char[] patientArr2) {

		for (int itr2 = 0; itr2 < patientArr1.length + 1; itr2++) {
			scoreMatrix[0][itr2] = 0;
		}
		System.out.println("Initialized...");
		Scoring(refArr, patientArr1, patientArr2);
	}

	public void Scoring(char[] refArr, char[] patientArr1, char[] patientArr2) {

		int scoreDiag = 0;
		int scoreUp = 0;
		int scoreLeft = 0;
		int curMaxScore = 0;
		int prevMaxScore = -1;
		int MaxX = 0;
		int MaxY = 0;

		int substitutionVal = 0;
		int counter = 1;// move in disease reference
		int handler = 0;// use to get up scores
		int mover = 0;// use to get left & right scores
		String path = "";

		while (refArr.length >= counter) {// X (through reference seq)

			if (counter % 2 == 0) {// even
				scoreMatrix[0][0] = 0;
				handler = 0;
				mover = 1;
			} else {// odd
				scoreMatrix[1][0] = 0;
				handler = 1;
				mover = 0;
			}
			path = "";// reset the path
			for (int itr2 = 1; itr2 <= patientArr1.length; itr2++) {// Y (through patient seq)

				substitutionVal = calculateSubstitution(refArr[counter - 1],
						patientArr1[itr2 - 1]);
				scoreDiag = scoreMatrix[mover][itr2 - 1] + substitutionVal;
				scoreUp = scoreMatrix[handler][itr2 - 1] + GAP_PANELTY;
				scoreLeft = scoreMatrix[mover][itr2] + GAP_PANELTY;
				curMaxScore = Math.max(scoreDiag,
						Math.max(scoreUp, scoreLeft));
				
				if(curMaxScore >= prevMaxScore){
					// store max score and position
					prevMaxScore = curMaxScore;
					MaxX = itr2-1;
					MaxY = counter-1;
				}
				
				scoreMatrix[handler][itr2] = calcScore(curMaxScore);
				//System.out.println("scoreMatrix PRINT"+scoreMatrix[handler][itr2]);
				path += storePath(handler, counter, itr2, scoreDiag, scoreUp,
						scoreLeft) + ",";
			}
			pathMatrix[counter - 1] = path;
			//System.out.println("PARTH PRINT"+path);
			//System.out.println("PARTH MaTrix"+pathMatrix[counter - 1]);
			counter++;
		}
		endIndex = MaxY;
		System.out.println("Scored...Completed");
		traceBack(refArr, patientArr1, patientArr2, MaxX, MaxY);
	}

	/**
	 * Calculate substitution value for Nucleotides
	 * @param ref
	 * @param target
	 * @return
	 */
	public int calculateSubstitution(char ref, char target) {

		if (ref == target) {

			return 1;
		} else {

			return -1;
		}
	}

	/**
	 * calculate the actual score
	 * @param scoreMax
	 * @return
	 */
	private int calcScore(int scoreMax){
		if(scoreMax < 0){
			return 0;
		}
		else{
			return scoreMax;
		}
	}
	
	/**
	 * Store the path to use in trace back
	 * @param handler
	 * @param counter
	 * @param itr2
	 * @param scoreDiag
	 * @param scoreUp
	 * @param scoreLeft
	 * @return
	 */
	public int storePath(int handler, int counter, int itr2, int scoreDiag,
			int scoreUp, int scoreLeft) {

		if (scoreMatrix[handler][itr2] == scoreDiag) {// store 0 (from diagonal)
			return 0;
		} else if (scoreMatrix[handler][itr2] == scoreUp) {// store 1 (from up)
			return 1;

		} else if (scoreMatrix[handler][itr2] == scoreLeft) {// store 2 (from left)
			return 2;

		} else {
			return -1;// nothing happens
		}

	}

	/**
	 * Trace back to get the alignment
	 * @param refArr
	 * @param patientArr1
	 * @param patientArr2
	 * @param MaxX
	 * @param MaxY
	 */
	public void traceBack(char[] refArr, char[] patientArr1, char[] patientArr2, int MaxX, int MaxY) {
		String ref = "";
		String target_1 = "";
		String target_2 = "";
		int nextX = MaxY;
		int nextY = MaxX;
		int refLen = MaxY;
		int pathIndex = nextY;

		int prevPath = Integer
				.parseInt(pathMatrix[nextX].split(",")[pathIndex]);

		while (nextX >= 0 && nextY >= 0) {

			if (prevPath == 0) {// diag
				ref += refArr[refLen];
				target_1 += patientArr1[nextY];//sequence 1
				target_2 += patientArr2[nextY];//sequence 2
				nextY = nextY - 1;
				refLen = refLen - 1;
				pathIndex = pathIndex - 1;
				nextX = nextX - 1;
				if (nextX >= 0 && pathIndex >= 0) {
					prevPath = Integer
							.parseInt(pathMatrix[nextX].split(",")[pathIndex]);
				}
			} else if (prevPath == 1) {// up
				ref += "_";
				target_1 += patientArr1[nextY];//sequence 1
				target_2 += patientArr2[nextY];//sequence 2
				nextY = nextY - 1;
				pathIndex = pathIndex - 1;

				if (nextX >= 0 && pathIndex >= 0) {
					prevPath = Integer
							.parseInt(pathMatrix[nextX].split(",")[pathIndex]);
				}
			} else if (prevPath == 2) {// left
				ref += refArr[refLen];
				target_1 += "_";//sequence 1
				target_2 += "_";//sequence 2
				refLen = refLen - 1;
				nextX = nextX - 1;

				if (nextX >= 0 && pathIndex >= 0) {
					prevPath = Integer
							.parseInt(pathMatrix[nextX].split(",")[pathIndex]);
				}
			} else {//score is -1
				break;
			}
			
			//System.out.println("pathMatrix[nextX]"+pathMatrix[nextX]);
		}

		ref1 = new StringBuffer(ref).reverse().toString();
		target1 = new StringBuffer(target_1).reverse().toString();
		target2 = new StringBuffer(target_2).reverse().toString();
		startIndex = nextX;
		System.out.println("Traceback completed...");

	}

	
	
	public void printPath(char[] refArr, char[] patientArr) {

		for (int itr2 = 0; itr2 < refArr.length; itr2++) {// Y
			System.out.print(pathMatrix[itr2] + ",");
			System.out.print("\n");
		}
	}


}


/*

import java.io.BufferedReader;
import java.io.FileReader;




public class GT_LocalAlignmentManager {

	private String refSeq = "";
//	String patient = "TCATACATTTTTCTCTAACTGCAAACATAATGTTTTCCCTTGTATTTTACAGATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTGGTAAAACCATTTGTTTTCTTCTTCTTCTTCTTCTTCTTTTCTTTTTTTTTTCTTTTTTTTTTTTGAGATGGAGTCTTGCTCTGTGGCCCAGGCTAGAAGCAGTCCTCCTGCCTTAGCCCCCTTAGTAGCTGGGATTACAGGCACGCGCCACCATGCCAGGCTAATTTTTGTATTTTTAGTAGAGACGGGGTTTCATCATGTTGGCCAGGCTGGTCTCGAACTCCTAACCTCAGGTGATCCACCCACCTCGGCTCCCCAAATTGCTGGGATTACAGGTGTGAGCCACTGTGCCCGGCCGGTAAAACCATTTTCATTTATTCTGGCAACATCTCTTTATGAGCATTGTGAATATGTTAGTGAATGTGCTAGATGCTCATAGATTTATATAAAAAGTTAGTGAAGAAGGAAAGATGGTATATTAAGTGGTTAGACAAGTGTTCTAATCAGTTAGAGTTCAGAGAAGGTCAGGGTACCTGATATAATCAAGAGAGAGACCTTACAGCCAGGTGAGGTGAATGTACCTATAATCCCAGCTACTTAGGAGGCTGAAATGGGAGGATCACTTGAGTCCAGGTTTGAGACCAGCCCAGGCAACATAGCAAGATCCCCATCAGATACACCAAAAAGACAGATTTCTTTTTTTTTTTTTTTTTTGAGACAGAGTCTCGCTCTGTCGCCCAGGCTGGAGCGCAGTGACACGATGTCAGCTCACTGCAACCTCCGCCTCCCAGGTTCAAGTGATTCTCCTGCCTCAGCCTCCTGAGTAGTTGGGACTACAGGGGTACGACACCAGACCTGGCTAATTTTTGTAATTTTAGTAGAGTCGGGGTTTCACCATATTGGTCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCCACCCTCCTTGGCCTCCCAGAGTGCTGGGATTACAGGCGTGAGCCACCAAGCCCGGCCAAAAAAGAGAGCTCTTATAGGCCCTTCCTTGCTTTGGAGCTTTATCTGCTCTGTGATGCTTATCTAAAATAGCCATAAGGTCACTGATATTTTTAAGCATTTGGAAATTACTTCAGCTGGGTGCCATGGCTCATGCCTATAATCCCAACCCTTTGGGAGGCTGA";
	private String patientSeq1 = "";
	private String patientSeq2 = "";

	private int[][] scoreMatrix = null;
	private String[] pathMatrix = null;
	private String ref1 = null;
	private String target1 = null;
	private String target2 = null;
	private int startIndex = 0;
	private int endIndex = 0;

	final int GAP_PANELTY = -1;
	
	public GT_LocalAlignmentManager (String patientSeq1, String patientSeq2,String diseaseId,int core,int fileNo){
		this.patientSeq1 = patientSeq1;
		this.patientSeq2 = patientSeq2;
		
		//start test data
		/*this.patientSeq1 = "TCATACATTTTTCTCTAACTGCAAACATAATGTTTTCCCTTGTATTTTACAGATGCAAAG"
			+ "AGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTAA"
			+ "CATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATC"
			+ "CTTCCTTGGTAAAACCATTTGTTTTCTTCTTCTTCTTCTTCTTCTTTTCTTTTTTTTTTCTT"
			+ "TTTTTTTTTTGAGATGGAGTCTTGCTCTGTGGCCCAGGCTAGAAGCAGTCCTCCTGCCTTAG"
			+ "CCCCCTTAGTAGCTGGGATTACAGGCACGCGCCACCATGCCAGGCTAATTTTTGTATTTTTA"
			+ "GTAGAGACGGGGTTTCATCATGTTGGCCAGGCTGGTCTCGAACTCCTAACCTCAGGTGATCC"
			+ "ACCCACCTCGGCTCCCCAAATTGCTGGGATTACAGGTGTGAGCCACTGTGCCCGGCCGGTAA"
			+ "AACCATTTTCATTTATTCTGGCAACATCTCTTTATTGAGCATTGTGAATATGTTAGTGAATG"
			+ "TGCTAGATGCTCATAGATTTATATAAAAAGTTAGTGAAGAAGGAAAGATGGTATATTAAGTG"
			+ "GTTAGACAAGTGTTCTAATCAGTTAGAGTTCAGAGAAGGTCAGGGTACCTGATATAATCAAG"
			+ "AGAGAGACCTTACAGCCAGGTGAGGTGAATGTACCTATAATCCCAGCTACTTAGGAGGCT";
		this.patientSeq2 = "TCATACATTTTTCTCTAACTGCAAACATAATGTTTTCCCTTGTATTTTACAGATGCAAAC"
		+ "AGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTAT"
		+ "CATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATC"
		+ "CTTCCTTGGTAAAACCATTTGTTTTCTTCTTCTTCTTCTTCTTCTTTTCTTTTTTTTTTCTT"
		+ "TTTTTTTTTTGAGATGGAGTCTTGCTCTGTGGCCCAGGCTAGAAGCAGTCCTCCTGCCTTAG"
		+ "CCCCCTTAGTAGCTGGGATTACAGGCACGCGCCACCATGCCAGGCTAATTTTTGTATTTTTA"
		+ "GTAGAGACGGGGTTTCATCATGTTGGCCAGGCTGGTCTCGAACTCCTAACCTCAGGTGATCC"
		+ "ACCCACCTCGGCTCCCCAAATTGCTGGGATTACAGGTGTGAGCCACTGTGCCCGGCCGGTAA"
		+ "AACCATTTTCATTTATTCTGGCAACATCTCTTTATTGAGCATTGTGAATATGTTAGTGAATG"
		+ "TGCTAGATGCTCATAGATTTATATAAAAAGTTAGTGAAGAAGGAAAGATGGTATATTAAGTG"
		+ "GTTAGACAAGTGTTCTAATCAGTTAGAGTTCAGAGAAGGTCAGGGTACCTGATATAATCAAG"
		+ "AGAGAGACCTTACAGCCAGGTGAGGTGAATGTACCTATAATCCCAGCTACTTAGGAGGCT";
		//end test data
		
		ref1 = null;
		target1 = null;
		target2 = null;
		startIndex = 0;
		this.LocalAlignment(diseaseId);
             
	}
	
	private void LocalAlignment(String diseaseId) {
		 this.refSeq = readReferenceFile(diseaseId);
		// patient = readReferenceFile();

		char[] refArr = this.refSeq.toCharArray();
		char[] patientArr1 = this.patientSeq1.toCharArray();
		char[] patientArr2 = this.patientSeq2.toCharArray();
		scoreMatrix = new int[2][patientArr1.length + 1];
		pathMatrix = new String[refArr.length];
		initialize(refArr, patientArr1, patientArr2);
//		printPath(refArr,patientArr);
	}

	private String readReferenceFile(String diseaseId) {
		String querySeq = null;
		try {
			BufferedReader br = new BufferedReader(new FileReader(diseaseId));
			String line = br.readLine();
			while (line != null) {
				if (line != null) {
					querySeq += line;
				}
				line = br.readLine();
			}
			br.close();
			System.out.println("File readed...");

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return querySeq;
	}

	public void initialize(char[] refArr, char[] patientArr1, char[] patientArr2) {

		for (int itr2 = 0; itr2 < patientArr1.length + 1; itr2++) {
			scoreMatrix[0][itr2] = 0;
		}

		Scoring(refArr, patientArr1, patientArr2);
	}

	public void Scoring(char[] refArr, char[] patientArr1, char[] patientArr2) {

		int scoreDiag = 0;
		int scoreUp = 0;
		int scoreLeft = 0;
		int curMaxScore = 0;
		int prevMaxScore = -1;
		int MaxX = 0;
		int MaxY = 0;

		int substitutionVal = 0;
		int counter = 1;
		int handler = 0;
		int mover = 0;
		String path = "";

		while (refArr.length >= counter) {// X
			System.out.println(100*counter/refArr.length+"%");
			if (counter % 2 == 0) {// even
				scoreMatrix[0][0] = 0;
				handler = 0;
				mover = 1;
			} else {// odd
				scoreMatrix[1][0] = 0;
				handler = 1;
				mover = 0;
			}
			path = "";// reset the path
			for (int itr2 = 1; itr2 <= patientArr1.length; itr2++) {// Y

				substitutionVal = calculateSubstitution(refArr[counter - 1],
						patientArr1[itr2 - 1]);
				scoreDiag = scoreMatrix[mover][itr2 - 1] + substitutionVal;
				scoreUp = scoreMatrix[handler][itr2 - 1] + GAP_PANELTY;
				scoreLeft = scoreMatrix[mover][itr2] + GAP_PANELTY;
				curMaxScore = Math.max(scoreDiag,
						Math.max(scoreUp, scoreLeft));
				
				if(curMaxScore >= prevMaxScore){
					//store max score and position
					prevMaxScore = curMaxScore;
					MaxX = itr2-1;
					MaxY = counter-1;
				}
				
				scoreMatrix[handler][itr2] = calcScore(curMaxScore);

				path += storePath(handler, counter, itr2, scoreDiag, scoreUp,
						scoreLeft) + ",";
				

			}
//			System.out.print("\n");
			pathMatrix[counter - 1] = path;
			counter++;
		}
//		System.out.println(MaxX+"::"+MaxY);
		endIndex = MaxY;
		traceBack(refArr, patientArr1, patientArr2, MaxX, MaxY);
	}

	public int calculateSubstitution(char ref, char target) {

		if (ref == target) {

			return 1;
		} else {

			return -1;
		}
	}

	private int calcScore(int scoreMax){
		if(scoreMax<0){
			return 0;
		}
		else{
			return scoreMax;
		}
	}
	
	public int storePath(int handler, int counter, int itr2, int scoreDiag,
			int scoreUp, int scoreLeft) {

		if (scoreMatrix[handler][itr2] == scoreDiag) {// store 0 (from diag)
			return 0;
		} else if (scoreMatrix[handler][itr2] == scoreUp) {// store 1 (from up)
			return 1;

		} else if (scoreMatrix[handler][itr2] == scoreLeft) {// store 2 (from
																// left)
			return 2;

		} else {
			return -1;// nothing happens
		}

	}

	public void traceBack(char[] refArr, char[] patientArr1, char[] patientArr2, int MaxX, int MaxY) {
		String ref = "";
		String target_1 = "";
		String target_2 = "";
		int nextX = MaxY;
		int nextY = MaxX;
		int refLen = MaxY;
		int pathIndex = nextY;

		int prevPath = Integer
				.parseInt(pathMatrix[nextX].split(",")[pathIndex]);

		while (nextX >= 0 && nextY >= 0) {

			if (prevPath == 0) {// diag
				ref += refArr[refLen];
				target_1 += patientArr1[nextY];
				target_2 += patientArr2[nextY];
				nextY = nextY - 1;
				refLen = refLen - 1;
				pathIndex = pathIndex - 1;
				nextX = nextX - 1;
				if (nextX >= 0 && pathIndex >= 0) {
					prevPath = Integer
							.parseInt(pathMatrix[nextX].split(",")[pathIndex]);
				}
			} else if (prevPath == 1) {// up
				ref += "_";
				target_1 += patientArr1[nextY];
				target_2 += patientArr2[nextY];
				nextY = nextY - 1;
				pathIndex = pathIndex - 1;

				if (nextX >= 0 && pathIndex >= 0) {
					prevPath = Integer
							.parseInt(pathMatrix[nextX].split(",")[pathIndex]);
				}
			} else if (prevPath == 2) {// left
				ref += refArr[refLen];
				target_1 += "_";
				target_2 += "_";
				refLen = refLen - 1;
				nextX = nextX - 1;

				if (nextX >= 0 && pathIndex >= 0) {
					prevPath = Integer
							.parseInt(pathMatrix[nextX].split(",")[pathIndex]);
				}
			} else {//score is 0
				break;
			}
		}
//		System.out.println(new StringBuffer(ref).reverse().toString());
//		System.out.println("************"+nextX);
		ref1 = new StringBuffer(ref).reverse().toString();
		target1 = new StringBuffer(target_1).reverse().toString();
		target2 = new StringBuffer(target_2).reverse().toString();
		startIndex = nextX;
	}

	public void printPath(char[] refArr, char[] patientArr) {

		for (int itr2 = 0; itr2 < refArr.length; itr2++) {// Y

			System.out.print(pathMatrix[itr2] + ",");

			System.out.print("\n");

		}

	}

	public String getRefSeq() {
		return refSeq;
	}

	public void setRefSeq(String refSeq) {
		this.refSeq = refSeq;
	}

	public String getPatientSeq1() {
		return patientSeq1;
	}

	public void setPatientSeq1(String patientSeq1) {
		this.patientSeq1 = patientSeq1;
	}

	public String getPatientSeq2() {
		return patientSeq2;
	}

	public void setPatientSeq2(String patientSeq2) {
		this.patientSeq2 = patientSeq2;
	}

	public String getRef1() {
		return ref1;
	}

	public void setRef1(String ref1) {
		this.ref1 = ref1;
	}

	public String getTarget1() {
		return target1;
	}

	public void setTarget1(String target1) {
		this.target1 = target1;
	}

	public String getTarget2() {
		return target2;
	}

	public void setTarget2(String target2) {
		this.target2 = target2;
	}


} */
