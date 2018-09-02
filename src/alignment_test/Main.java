
/**
 *  Local Alignment with Bio Java Package
 *  @Genotec
 */

package alignment_test;

import java.io.File;


public class Main {

    /**
     * @param args the command line arguments
     */
    static double d = System.currentTimeMillis();
    static int count = 4;
    static int perCore = 1;

    public static void main(String[] args) throws Exception {
    	 int me = 0;
        try {
        	 // Dummy patient sequence taken from BRACA Gene (1)
        	String dummyPatientBRACASeq1 = "TCATACATTTTTCTCTAACTGCAAACATAATGTTTTCCCTTGTATTTTACAGATGCAAAG"
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
        	
    		
    	//String referanceFilePath = "C:/Users/Win 8/Desktop/final year project re write/BRCA1_REF.txt";
        	// Path of the Refernace file BRCA
        	String referanceFilePath = "<Path of File>/BRCA1_REF.txt";
        	SequenceAlignmentManager.alignSequenceLocal(dummyPatientBRACASeq1,referanceFilePath);
        } catch (Exception ex) {
            System.out.print(ex);
        }
        System.out.println("Core No:" + me + " done____________________");
        System.out.println(System.currentTimeMillis() - d);

    }
    private static final long MEGABYTE = 1024L * 1024L;

    public static long bytesToMegabytes(long bytes) {
        return bytes / MEGABYTE;
    }
}
