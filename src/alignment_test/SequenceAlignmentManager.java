/**
 *  Local Alignment with Bio Java Package
 *  @Genotec
 */
package alignment_test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;

//import com.genotec.util.Constants;

public class SequenceAlignmentManager {
public static SequencePair<DNASequence, NucleotideCompound> alignSequenceLocal(String sequence,String refPath){
		
		String targetSeq = sequence;
		String querySeq="";


		File file = new File(refPath);
		try {
			BufferedReader br= new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while (line != null)
			{
				if(line!=null){
				querySeq +=line;
				}
				line = br.readLine();
			   System.out.println(querySeq);
		
			}
			br.close();
			System.out.println("File readed...");
		
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		DNASequence query;
		try {
			System.out.println("Alignment Starting!");
			DNASequence target;
			target = new DNASequence(targetSeq,
					AmbiguityDNACompoundSet.getDNACompoundSet());
			
			query = new DNASequence(querySeq,
					AmbiguityDNACompoundSet.getDNACompoundSet());
			SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

			SimpleGapPenalty gapP = new SimpleGapPenalty();
			gapP.setOpenPenalty((short)5);
			gapP.setExtensionPenalty((short)2);

			SequencePair<DNASequence, NucleotideCompound> psa =
					Alignments.getPairwiseAlignment(query, target,
							PairwiseSequenceAlignerType.LOCAL, gapP, matrix);

			PairwiseSequenceAligner<DNASequence, NucleotideCompound> psa2 =
				Alignments.getPairwiseAligner(query, target,
						PairwiseSequenceAlignerType.LOCAL, gapP, matrix);

			System.out.println(psa);
			return psa;

		} catch (CompoundNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}

		
	}
	/*public static SequencePair<DNASequence, NucleotideCompound> alignSequenceGlobal(String sequence){
		
		String targetSeq = sequence;
		DNASequence target;
		try {
			target = new DNASequence(targetSeq,
					AmbiguityDNACompoundSet.getDNACompoundSet());
		} catch (CompoundNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		String querySeq="";
		File file = new File("BRCA1_REF.txt");
		try {
			BufferedReader br= new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while (line != null)
			{
				if(line!=null){
				querySeq +=line;
				}
				line = br.readLine();		
			}
			br.close();
			System.out.println("File readed...");
		
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		DNASequence query;
		try {
			query = new DNASequence(querySeq,
					AmbiguityDNACompoundSet.getDNACompoundSet());
		} catch (CompoundNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

		SimpleGapPenalty gapP = new SimpleGapPenalty();
		gapP.setOpenPenalty((short)5);
		gapP.setExtensionPenalty((short)2);

		SequencePair<DNASequence, NucleotideCompound> psa =
				Alignments.getPairwiseAlignment(query, target,
						PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);

		System.out.println(psa);
		return psa;

	} */
}
