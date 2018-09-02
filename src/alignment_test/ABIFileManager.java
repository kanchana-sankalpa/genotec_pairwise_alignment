package alignment_test;

import java.io.File;


public class ABIFileManager {


	private char hetroSequence[];
	private String sequence1, sequence2 = "";


	public void getSequences(File file,double hetroLimit) throws Exception {

		sequence1 = "";
		sequence2 = "";
		ABITrace trace = new ABITrace(file,hetroLimit);

		hetroSequence = trace.getHetroSequence();
		sequence1 = trace.getSequncetrace();
		for (int i = 0; i < sequence1.length(); i++) {
			if (hetroSequence[i] != '\u0000') {
				sequence2 = sequence2 + hetroSequence[i];
			} else {
				sequence2 = sequence2 + sequence1.charAt(i);
			}
		}

	}


	public String getSequence1() {
		return sequence1;
	}


	public void setSequence1(String sequence1) {
		this.sequence1 = sequence1;
	}


	public String getSequence2() {
		return sequence2;
	}


	public void setSequence2(String sequence2) {
		this.sequence2 = sequence2;
	}

}