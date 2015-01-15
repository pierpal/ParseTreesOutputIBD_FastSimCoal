import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;


/**
 * Parses Newick trees and outputs IBD segments
 * @author Pier Palamara <pier@cs.columbia.edu>
 */

public class ParseTreesOutputIBD_FastSimCoal {

    static double factor;
    static int end;
    static InputStreamReader stream = new InputStreamReader(System.in);
    static final boolean germline = false;
    static int trees = 0;
    static boolean outputProgress = false;
    static final int howOften = 5000;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        //	args[0] is the threshold line, args[1] is the blocksize, args[2] is the minLength

        factor = Double.valueOf(args[0]);
        end = Integer.valueOf(args[1]);
        if (args.length == 4) {
            try {
                int t = Integer.parseInt(args[3]);
                trees = t;
                outputProgress = true;
            } catch (NumberFormatException E) {
                FileInputStream fstream = new FileInputStream(args[3]);
                DataInputStream in = new DataInputStream(fstream);
                stream = new InputStreamReader(in);
                System.err.println("Reading trees from file " + args[3]);
            }
        }

        System.err.println("Mb/cM: " + factor + "\tend: " + end + "\tminimum length: " + Float.parseFloat(args[2]));
        getSegments(Float.parseFloat(args[2]));

    }

    static void getSegments(float minLen) throws IOException {
        //Takes in a file and puts it into a buffer.
        int position = 0;
        //This Arraylist of an integer array stores each array of integers all sharing SNP.
        ArrayList<Integer[]> allStartsOverall = new ArrayList<Integer[]>();
        //This string is used to read in the input

        /*This reads in the input,
         * parses each individual block tree,
         * inputs the string from each tree into c1.getSets() method
         * inputs an arraylist of sets of IBD integers called list to getSets() to add, an arraylist of SNP individuals
         * for each set of IBD integers, getSets will perform the method out()
         *
         */

        String str = "";
        BufferedReader in = new BufferedReader(stream);
        str = in.readLine();


        String[][] ancestor = null;
        String[][] past_ancestor = null;
        Integer[][] IBD = null;

        int counter = 0;
        int N = 0;

        long start = System.currentTimeMillis(), time = start;

        while (true) {
            str = in.readLine();
            if (str == null) {
                break;
            }
            String[] splitString = str.split("\\s+");
            if (splitString.length > 1 && !splitString[1].isEmpty()) {
                counter++;
                if (counter % howOften == 0) {
                    time = System.currentTimeMillis();
                    double frac = counter/(double)trees;
                    double ETA=(time-start)/frac*(1-frac);
                    if (outputProgress) {
                        System.err.format("Progress:\t%.2f\tETA: \t%d:%d:%d\n", (100.0*frac), (int)ETA/1000/60/60, ((int)ETA/1000/60)%60, ((int)ETA/1000)%60);
                    }
                }
                position = Integer.valueOf(splitString[0]) + 1;
//                System.out.println(position / 1000000.0);
                String s = "Begin Trees;\nTree Tree0= \n" + splitString[1] + "\nEnd";
                ArrayList<ArrayList<Integer>> list = new ArrayList<ArrayList<Integer>>();
                Sets c1 = new Sets();
                ancestor = c1.getSets(s);
                if (counter == 1) {
                    N = Sets.id2int.size();
                    IBD = new Integer[N][N];
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < N; j++) {
                            IBD[i][j] = 0;
                        }
                    }
                    past_ancestor = ancestor.clone();
                }
//                System.out.println("N is " + N);
                for (int i = 0; i < N; i++) {
                    for (int j = i + 1; j < N; j++) {
                        //System.out.println("CG " + ancestor[i][j].split("\\.")[0]);
                        if (!ancestor[i][j].equalsIgnoreCase(past_ancestor[i][j])) {
//                                System.out.println("ancestor of " + i + " " + j + " has changed from " + ancestor[i][j] + " to " + past_ancestor[i][j]);
                            double len = Math.round((((double) (position - IBD[i][j])) * factor)) / 1000000.0;
                            if (len >= minLen) {
                                if (!germline) {
//                                    System.out.println(Sets.int2id.get(i) + "\t" + Sets.int2id.get(j) + "\t" + (IBD[i][j] + 1) + "\t" + position + "\t" + len + "\t" + past_ancestor[i][j]);
                                    System.out.println(Sets.int2id.get(i) + "\t" + Sets.int2id.get(j) + "\t" + (IBD[i][j] + 1) + "\t" + position + "\t" + len );
                                } else {
                                    System.out.println(Sets.int2id.get(i) + "\t" + Sets.int2id.get(i) + "\t" + Sets.int2id.get(j) + "\t" + Sets.int2id.get(j) + "\tNA\t" + (IBD[i][j] + 1) + "\t" + position + "\tNA\tNA\tNA\t" + len + "\tcM\t" + past_ancestor[i][j] + "\t" + 0 + "\t" + 0);
                                }
                            }
                            IBD[i][j] = position;
                        } // end if
                    } // end for j
                } // end for i
                past_ancestor = ancestor.clone();
            }
        }

        position = end;
        // at end of chromosome print all
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double len = Math.round((((double) (position - IBD[i][j])) * factor)) / 1000000.0;
                if (len >= minLen) {
                    if (!germline) {
//                        System.out.println(Sets.int2id.get(i) + "\t" + Sets.int2id.get(j) + "\t" + (IBD[i][j] + 1) + "\t" + position + "\t" + len + "\t" + past_ancestor[i][j]);
                        System.out.println(Sets.int2id.get(i) + "\t" + Sets.int2id.get(j) + "\t" + (IBD[i][j] + 1) + "\t" + position + "\t" + len);
                    } else {
                        System.out.println(Sets.int2id.get(i) + "\t" + Sets.int2id.get(i) + "\t" + Sets.int2id.get(j) + "\t" + Sets.int2id.get(j) + "\tNA\t" + (IBD[i][j] + 1) + "\t" + position + "\tNA\tNA\tNA\t" + len + "\tcM\t" + past_ancestor[i][j] + "\t" + 0 + "\t" + 0);
                    }
                }
            }
        }
    }
}
