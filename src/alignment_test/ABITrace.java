/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package alignment_test;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Title: ABITrace<br>
 * <br>
 * ABITrace is a class for managing ABI file information, it is capable of
 * opening an ABI file and storing the most important fields, which can be
 * recalled as simple java types. It can also return an image corresponding to
 * the trace. It has three constructors with input types
 * <code>File, URL, and byte[]</code>.<br>
 * <br>
 * ABI files contain two sets of base call and sequence data, one that was
 * originally created programmatically and the other, which is an editable copy.
 * This version of this object only references the original unedited data.<br>
 * 
 * Copyright (c) 2001
 * 
 * @author David H. Klatte, Ph.D.
 * @author Matthew Pocock
 * @version 0.5alpha Edited by Waruna Jayaweera
 */
public class ABITrace {

	private String sequence, sequncetrace = "";
	private int A[], G[], C[], T[], Basecalls[];
	private char hetroChannels[], hetroSequence[], sequnceChannel[];
	private HashMap<Integer, Integer> extremeOutliers = new HashMap<Integer,Integer>();
	private int TraceLength, SeqLength;
	private byte[] TraceData;
	private int MacJunk = 0, IndexBase, PLOC;
	private  int AbsIndexBase = 26; // The file location of the Index
											// pointer
	// the next declaration is for the actual file pointers
	private int DATA9, DATA10, DATA11, DATA12, PBAS2, FWO;
	private int Q1, Q2, Q3,hetroCount,minimum,maximum,IQR , OFL , OFR ;

	public ABITrace(File ABIFile,double hetroLimit) throws IOException {

		byte[] bytes = null;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		FileInputStream fis = new FileInputStream(ABIFile);
		BufferedInputStream bis = new BufferedInputStream(fis);
		int b;
		while ((b = bis.read()) >= 0) {
			baos.write(b);
		}
		bis.close();
		fis.close();
		baos.close();
		bytes = baos.toByteArray();
		initData(bytes, hetroLimit);
	}

	/**
	 * The URL constructor opens an ABI file from any URL.
	 * 
	 * @param ABIFile
	 *            is a <code>java.net.URL</code> for an ABI trace file.
	 * @throws IOException
	 *             if there is a problem reading from the URL.
	 * @throws IllegalArgumentException
	 *             if the URL does not contain a valid ABI file.
	 */
	public ABITrace(URL ABIFile,double hetroLimit) throws IOException {
		byte[] bytes = null;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		InputStream is = ABIFile.openStream();
		BufferedInputStream bis = new BufferedInputStream(is);
		int b;
		while ((b = bis.read()) >= 0) {
			baos.write(b);
		}
		bis.close();
		is.close();
		baos.close();
		bytes = baos.toByteArray();
		initData(bytes,hetroLimit);
	}

	/**
	 * The <code>byte[]</code> constructor parses an ABI file represented as a
	 * byte array.
	 * 
	 * @throws IllegalArgumentException
	 *             if the data does not represent a valid ABI file.
	 */
	public ABITrace(byte[] ABIFileData,int hetroLimit) {
		initData(ABIFileData,hetroLimit);
	}

	/**
	 * Returns the length of the sequence (number of bases) in this trace.
	 */
	public int getSequenceLength() {
		return SeqLength;
	}

	/**
	 * Returns the length of the trace (number of x-coordinate points in the
	 * graph).
	 */
	public int getTraceLength() {
		return TraceLength;
	}

	/**
	 * Returns an <code>int[]</code> array that represents the base calls - each
	 * int in the array corresponds to an x-coordinate point in the graph that
	 * is a peak (a base location).
	 */
	public int[] getBasecalls() {
		return Basecalls;
	}

	private void initData(byte[] fileData,double hetroLimit) {
		TraceData = fileData;
		if (isABI()) {
			setIndex();
			setBasecalls();
			setSeq();
			setTraces();
			sethetroChannelPositions(hetroLimit);
		} else {
			throw new IllegalArgumentException("Not a valid ABI file.");
		}
	}

	/**
	 * A utility method which fills array b with data from the trace starting at
	 * traceDataOffset.
	 */
	private void getSubArray(byte[] b, int traceDataOffset) {
		for (int x = 0; x <= b.length - 1; x++) {
			b[x] = TraceData[traceDataOffset + x];
		}
	}

	/**
	 * Shuffle the pointers to point to the proper spots in the trace, then load
	 * the traces into their arrays.
	 */
	public void sethetroChannelPositions(double hetroLimit) {

		
		int[] peakValues = new int[SeqLength];
		hetroChannels = new char[TraceLength];
		sequnceChannel = new char[TraceLength];
		hetroSequence = new char[SeqLength];
		int peak = 0;
		char[] seq = sequence.toCharArray();
		char s = '\u0000';
		char c = ' ';
		for (int i = 0; i < SeqLength; i++) {
			c = seq[i];
			switch (c) {
			case 'A':
				sequncetrace = sequncetrace + c;
				peakValues[i] = A[Basecalls[i]];
				break;
			case 'C':
				sequncetrace = sequncetrace + c;
				peakValues[i] = C[Basecalls[i]];
				break;
			case 'G':
				sequncetrace = sequncetrace + c;
				peakValues[i] = G[Basecalls[i]];
				break;
			case 'T':
				sequncetrace = sequncetrace + c;
				peakValues[i] = T[Basecalls[i]];

				break;

			}

			if (c != 'A' && A[Basecalls[i]] < peakValues[i] && A[Basecalls[i]] > peakValues[i]*hetroLimit
					&& A[Basecalls[i] - 1] < A[Basecalls[i]]
					&& A[Basecalls[i]] > A[Basecalls[i] + 1]) {
				s = 'A';
				peak = A[Basecalls[i]];
			}
			if (c != 'C' &&C[Basecalls[i]] < peakValues[i] && C[Basecalls[i]] > peakValues[i]*hetroLimit
					&& C[Basecalls[i] - 1] < C[Basecalls[i]]
					&& C[Basecalls[i]] > C[Basecalls[i] + 1]) {
				if (C[Basecalls[i]] > peak) {
					peak = C[Basecalls[i]];
					s = 'C';
				}

			}
			if (c != 'G' &&G[Basecalls[i]] < peakValues[i] && G[Basecalls[i]] > peakValues[i]*hetroLimit
					&& G[Basecalls[i] - 1] < G[Basecalls[i]]
					&& G[Basecalls[i]] > G[Basecalls[i] + 1]) {
				if (G[Basecalls[i]] > peak) {
					peak = G[Basecalls[i]];
					s = 'G';
				}
			}
			if (c != 'T' &&T[Basecalls[i]] < peakValues[i] && T[Basecalls[i]] > peakValues[i]*hetroLimit
					&& T[Basecalls[i] - 1] < T[Basecalls[i]]
					&& T[Basecalls[i]] > T[Basecalls[i] + 1]) {
				if (T[Basecalls[i]] > peak) {
					s = 'T';
				}
			}
			if (s != '\u0000') {

				hetroChannels[Basecalls[i]] = s;
				hetroSequence[i] = s;
				hetroCount++;
				s = '\u0000';
				peak = 0;
			}
			sequnceChannel[Basecalls[i]] = c;

		}
		Arrays.sort(peakValues);
		minimum=peakValues[0];
		maximum=peakValues[peakValues.length-1];
		Q1 = (peakValues[SeqLength / 4] + peakValues[(SeqLength / 4) + 1]) / 2;
		Q2 = (peakValues[SeqLength / 2] + peakValues[(SeqLength / 2) + 1]) / 2;
		Q3 = (peakValues[SeqLength / 4 * 3] + peakValues[(SeqLength / 4 * 3) + 1]) / 2;
		IQR = Q3 - Q1;
		OFL = (int) (Q1 - (IQR * 3));
		OFR = (int) (Q3 + (IQR * 3));


		for (int i = 0; i < SeqLength; i++) {
			c = seq[i];
			switch (c) {
			case 'A':
				if (A[Basecalls[i]] > OFR || A[Basecalls[i]] < OFL) {
					extremeOutliers.put(i, A[Basecalls[i]]);
				}
				break;
			case 'C':
				if (C[Basecalls[i]] > OFR || C[Basecalls[i]] < OFL) {
					extremeOutliers.put(i, C[Basecalls[i]]);
				}
				break;
			case 'G':
				if (G[Basecalls[i]] > OFR || G[Basecalls[i]] < OFL) {
					extremeOutliers.put(i, G[Basecalls[i]]);
				}
				break;
			case 'T':
				if (T[Basecalls[i]] > OFR || T[Basecalls[i]] < OFL) {
					extremeOutliers.put(i, T[Basecalls[i]]);
				}
				break;

			}
		}
	}

	private void setTraces() {
		int pointers[] = new int[4]; // alphabetical, 0=A, 1=C, 2=G, 3=T
		int datas[] = new int[4];
		char order[] = new char[4];

		datas[0] = DATA9;
		datas[1] = DATA10;
		datas[2] = DATA11;
		datas[3] = DATA12;

		for (int i = 0; i <= 3; i++) {
			order[i] = (char) TraceData[FWO + i];
		}
		for (int i = 0; i <= 3; i++) {
			switch (order[i]) {
			case 'A':
			case 'a':
				pointers[0] = datas[i];
				break;
			case 'C':
			case 'c':
				pointers[1] = datas[i];
				break;
			case 'G':
			case 'g':
				pointers[2] = datas[i];
				break;
			case 'T':
			case 't':
				pointers[3] = datas[i];
				break;
			default:
				throw new IllegalArgumentException(
						"Trace contains illegal values.");
			}
		}
		A = new int[TraceLength];
		C = new int[TraceLength];
		G = new int[TraceLength];
		T = new int[TraceLength];

		for (int i = 0; i <= 3; i++) {
			byte[] qq = new byte[TraceLength * 2];
			getSubArray(qq, pointers[i]);
			DataInputStream dis = new DataInputStream(new ByteArrayInputStream(
					qq));
			for (int x = 0; x <= TraceLength - 1; x++) {
				try {
					if (i == 0) {
						A[x] = (int) dis.readShort();

					}
					if (i == 1) {
						C[x] = (int) dis.readShort();
					}
					if (i == 2) {
						G[x] = (int) dis.readShort();
					}
					if (i == 3) {
						T[x] = (int) dis.readShort();
					}

				} catch (IOException e)// This shouldn't happen. If it does
										// something must be seriously wrong.
				{
					throw new IllegalStateException(
							"Unexpected IOException encountered while manipulating internal streams.");
				}
			}

		}
	}

	/**
	 * Fetch the sequence from the trace data.
	 */
	private void setSeq() {
		char tempseq[] = new char[SeqLength];
		for (int x = 0; x <= SeqLength - 1; ++x) {
			tempseq[x] = (char) TraceData[PBAS2 + x];
		}
		sequence = new String(tempseq);
	}

	/**
	 * Fetch the base calls from the trace data.
	 */
	private void setBasecalls() {
		Basecalls = new int[SeqLength];
		byte[] qq = new byte[SeqLength * 2];
		getSubArray(qq, PLOC);
		DataInputStream dis = new DataInputStream(new ByteArrayInputStream(qq));
		for (int i = 0; i <= SeqLength - 1; ++i) {
			try {
				int a = (int) dis.readShort();
				;
				Basecalls[i] = a;
				// System.out.println(a+"|");
			} catch (IOException e)// This shouldn't happen. If it does
									// something must be seriously wrong.
			{
				throw new IllegalStateException(
						"Unexpected IOException encountered while manipulating internal streams.");
			}
		}
	}

	/**
	 * Utility method to return an int beginning at <code>pointer</code> in the
	 * TraceData array.
	 */
	private int getIntAt(int pointer) {
		int out = 0;
		byte[] temp = new byte[4];
		getSubArray(temp, pointer);
		try {
			DataInputStream dis = new DataInputStream(new ByteArrayInputStream(
					temp));
			out = dis.readInt();
		} catch (IOException e) // This shouldn't happen. If it does something
								// must be seriously wrong.
		{
			throw new IllegalStateException(
					"Unexpected IOException encountered while manipulating internal streams.");
		}
		return out;
	}

	/**
	 * Utility method to translate y coordinates from graph space (where up is
	 * greater) to image space (where down is greater).
	 */


	/**
	 * Get the maximum height of any of the traces. The data is persisted for
	 * performance in the event of multiple calls, but it initialize lazily.
	 */



	// calculates the necessary scaling to allow the trace to fit vertically
	// in the space specified.
	/**
	 * Returns the scaling factor necessary to allow all of the traces to fit
	 * vertically into the specified space.
	 * 
	 * @param <code>height</code> - the required height in pixels.
	 */


	/**
	 * Sets up all of the initial pointers to the important records in
	 * TraceData.
	 */
	private void setIndex() {
		int DataCounter, PBASCounter, PLOCCounter, NumRecords;
		byte[] RecNameArray = new byte[4];
		String RecName;

		DataCounter = 0;
		PBASCounter = 0;
		PLOCCounter = 0;

		IndexBase = getIntAt(AbsIndexBase + MacJunk);
		// System.out.println("IndexBase " + IndexBase);
		NumRecords = getIntAt(AbsIndexBase - 8 + MacJunk);
		// System.out.println("NumRecords " + NumRecords);

		for (int record = 0; record <= NumRecords - 1; record++) {
			getSubArray(RecNameArray, (IndexBase + (record * 28)));
			RecName = new String(RecNameArray);
			// System.out.println(RecName);
			if (RecName.equals("FWO_")) {
				FWO = IndexBase + (record * 28) + 20;
			}
			if (RecName.equals("DATA")) {
				++DataCounter;
				if (DataCounter == 9) {
					DATA9 = IndexBase + (record * 28) + 20;
				}
				if (DataCounter == 10) {
					DATA10 = IndexBase + (record * 28) + 20;
				}
				if (DataCounter == 11) {
					DATA11 = IndexBase + (record * 28) + 20;
				}
				if (DataCounter == 12) {
					DATA12 = IndexBase + (record * 28) + 20;
				}
			}
			if (RecName.equals("PBAS")) {
				++PBASCounter;
				if (PBASCounter == 2) {
					PBAS2 = IndexBase + (record * 28) + 20;
				}
			}
			if (RecName.equals("PLOC")) {
				++PLOCCounter;
				if (PLOCCounter == 2) {
					PLOC = IndexBase + (record * 28) + 20;
				}
			}

		} // next record
		TraceLength = getIntAt(DATA12 - 8);
		SeqLength = getIntAt(PBAS2 - 4);
		PLOC = getIntAt(PLOC) + MacJunk;
		DATA9 = getIntAt(DATA9) + MacJunk;
		DATA10 = getIntAt(DATA10) + MacJunk;
		DATA11 = getIntAt(DATA11) + MacJunk;
		DATA12 = getIntAt(DATA12) + MacJunk;
		PBAS2 = getIntAt(PBAS2) + MacJunk;
		// System.out.println(TraceLength + " " + SeqLength + " " + PLOC + " " +
		// DATA9 + " " + DATA10 + " " + DATA11 + " " + DATA12 + " " + PBAS2);
	}

	/**
	 * Test to see if the file is ABI format by checking to see that the first
	 * three bytes are "ABI". Also handle the special case where 128 bytes were
	 * prepended to the file due to binary FTP from an older macintosh system.
	 */
	private boolean isABI() {
		char ABI[] = new char[4];

		for (int i = 0; i <= 2; i++) {
			ABI[i] = (char) TraceData[i];
		}
		if (ABI[0] == 'A' && (ABI[1] == 'B' && ABI[2] == 'I')) {
			return true;
		} else {
			for (int i = 128; i <= 130; i++) {
				ABI[i] = (char) TraceData[i];
			}
			if (ABI[0] == 'A' && (ABI[1] == 'B' && ABI[2] == 'I')) {
				MacJunk = 128;
				return true;
			} else {
				return false;
			}
		}
	}

	public String getSequncetrace() {
		return sequncetrace;
	}

	public int[] getA() {
		return A;
	}

	public int[] getC() {
		return C;
	}

	public int[] getG() {
		return G;
	}

	public int[] getT() {
		return T;
	}

	public char[] getHetroChannels() {
		return hetroChannels;
	}

	public char[] getHetroSequence() {
		return hetroSequence;
	}

	public char[] getSequnceChannel() {
		return sequnceChannel;
	}

	public HashMap<Integer, Integer> getextremeOutliers() {
		return extremeOutliers;
	}

	public int getQ1() {
		return Q1;
	}

	public int getQ2() {
		return Q2;
	}

	public int getQ3() {
		return Q3;
	}

	public int getHetroCount() {
		return hetroCount;
	}

	public void setHetroCount(int hetroCount) {
		this.hetroCount = hetroCount;
	}

	public int getMinimum() {
		return minimum;
	}

	public void setMinimum(int minimum) {
		this.minimum = minimum;
	}

	public int getMaximum() {
		return maximum;
	}

	public void setMaximum(int maximum) {
		this.maximum = maximum;
	}

	public int getOFL() {
		return OFL;
	}

	public void setOFL(int oFL) {
		OFL = oFL;
	}

	public int getOFR() {
		return OFR;
	}

	public void setOFR(int oFR) {
		OFR = oFR;
	}
	

}
