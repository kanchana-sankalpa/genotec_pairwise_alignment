package alignment_test;

import java.io.File;

public class SequncetialAligner {
public static ABIFileManager am = new ABIFileManager();
    public static void main(String args[]) throws Exception {
        double d = System.currentTimeMillis();

            for (int j =1; j <= 2; j++) {

                am = new ABIFileManager();
                am.getSequences(new File("G:/" + j + ".ab1"), 0.25);
               /* GT_LocalAlignmentManager alignmentManager = new GT_LocalAlignmentManager(
                        am.getSequence1(), am.getSequence2(), "G:/BRCA1_REF.txt", 1, j); */
                System.out.println(j / 8 + "% done...______________");
         

            

        }
         
        System.out.println(System.currentTimeMillis() - d);
    }
    private static final long MEGABYTE = 1024L * 1024L;

    public static long bytesToMegabytes(long bytes) {
        return bytes / MEGABYTE;
    }
}
