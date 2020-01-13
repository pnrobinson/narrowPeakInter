package org.jax.npi.analysis;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jax.npi.data.RegulatoryElement;
import org.jax.npi.data.H3K27AcSignal;

import static org.apache.commons.math3.stat.inference.TestUtils.chiSquare;
import static org.apache.commons.math3.stat.inference.TestUtils.chiSquareTest;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class ChromosomeWithEnhancers {

    private final Map<String, List<RegulatoryElement>> chromosome2enhancerList;

    private int number_of_experiments;

    public ChromosomeWithEnhancers(List<RegulatoryElement> enhancerList) {
        chromosome2enhancerList = new HashMap<>();
        number_of_experiments = 0;
        for (RegulatoryElement e : enhancerList) {
            chromosome2enhancerList.putIfAbsent(e.getChromosome(), new ArrayList<>());
            List<RegulatoryElement> enlst = chromosome2enhancerList.get(e.getChromosome());
            enlst.add(e);
        }
        // sort the Enhancers
        for (List<RegulatoryElement> lst : chromosome2enhancerList.values()) {
            Collections.sort(lst);
        }
    }


    /**
     *     [0] string chrom;        "Reference sequence chromosome or scaffold"
     *     [1] uint   chromStart;   "Start position in chromosome"
     *     [2] uint   chromEnd;     "End position in chromosome"
     *     [3] string name;
     *     [4] uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
     *     [5] char[1]  strand;     "+ or - or . for unknown"
     *     [6] float  signalValue;  "Measurement of average enrichment for the region"
     *     [7] float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
     *     [8] float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
     *     [9] int   peak;         "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."
     * @param bedfile
     */
    public void addDataFromBedFile(File bedfile){
        try {
            InputStream fileStream = new FileInputStream(bedfile);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream);
            BufferedReader br = new BufferedReader(decoder);
            String line;
            while ((line = br.readLine())!=null) {
                String []F = line.split("\t");
                if (F.length != 10) {
                    throw new RuntimeException("Malformed BED file line: " + line);
                }
                String chrom = F[0];
                int begin = Integer.parseInt(F[1]);
                int end = Integer.parseInt(F[2]);
                double value = Double.parseDouble(F[6]);
                addDataPoint(chrom, begin, end, value);
            }
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("Could not read " + bedfile.getAbsolutePath());
        }
        this.number_of_experiments += 1;
        // process the data from this experment
        for ( List<RegulatoryElement> lst : chromosome2enhancerList.values()) {
            for (RegulatoryElement e : lst) {
                e.processLastExperiment();
            }
        }
    }



    /** Naive algorithm, but fast enough in practice. */
    private void addDataPoint(String chrom, int beginH3K27ac, int endH3K27ac, double value) {
        if (! this.chromosome2enhancerList.containsKey(chrom)) {
            if (chrom.contains("random") || chrom.contains("Un_")) {
                // do not worry about these scaffolds, just skip
                return;
            }
            throw new RuntimeException("Could not find chromosome " + chrom);
        }
        if (beginH3K27ac > endH3K27ac) {
            throw new RuntimeException(String.format("Begin=%d and end =%d\n", beginH3K27ac, endH3K27ac));
        }
        List<RegulatoryElement> enhancers = this.chromosome2enhancerList.get(chrom);
        for (RegulatoryElement e : enhancers) {
            int enhancerStartPos = e.getBegin();
            int enhancerEndPos = e.getEnd();
            if (enhancerStartPos > endH3K27ac) {
                // the enhancers are sorted -- the current enhancer is already beyond the end
                // of the new interval and there was no match. There is no more hope
                // of finding a match, so we can just return
                return;
            } else if (enhancerEndPos < beginH3K27ac) {
                // The H3K27ac region is located 3' (downstream) to the current enhancer
                // therefore, skip to the next enhancer and check for overlap
                continue;
            } else if (beginH3K27ac <= enhancerStartPos && endH3K27ac >= enhancerStartPos) {
                /*
                        |--enhancer--|
                     |------ peak---------|
                   begin                 end
                        or
                          |--enhancer--------|
                     |------ peak---------|
                   begin                 end
                 */
                int B = enhancerStartPos; //
                int E = Math.min(endH3K27ac, enhancerEndPos);
                H3K27AcSignal h3k27 = new H3K27AcSignal(B, E, value);
                e.addH3K27AcValue(h3k27);
            } else if (beginH3K27ac <= enhancerEndPos && endH3K27ac >= enhancerEndPos) {
                  /*
                        |--enhancer--|
                     |------ peak---------|
                   begin                 end
                        or
                 |--enhancer--------|
                     |------ peak---------|
                   begin                 end
                 */
                int B = beginH3K27ac < enhancerStartPos ? enhancerStartPos : beginH3K27ac;
                int E = endH3K27ac > enhancerEndPos ? enhancerEndPos : endH3K27ac;
                H3K27AcSignal h3k27 = new H3K27AcSignal(B, E, value);
                e.addH3K27AcValue(h3k27);
            } else if (enhancerStartPos <= beginH3K27ac &&  endH3K27ac <= enhancerEndPos) {
                // H3K27Ac region including entirely within the enhancer
                 /*
                        |----enhancer----|
                             |-- peak-|
                            begin     end
                */
                int B = Math.max(enhancerStartPos, beginH3K27ac);
                int E = Math.min(enhancerEndPos, endH3K27ac);
                H3K27AcSignal h3k27 = new H3K27AcSignal(B, E, value);
                e.addH3K27AcValue(h3k27);
            } else if ( beginH3K27ac <=  enhancerStartPos &&  enhancerEndPos  <= endH3K27ac) {
                // enhancer  including entirely within the H3K27Ac region
                 /*
                        |----enhancer----|
                      |------- peak----------------|
                */
                int B = Math.max(enhancerStartPos, beginH3K27ac);
                int E = Math.min(enhancerEndPos, endH3K27ac);
                H3K27AcSignal h3k27 = new H3K27AcSignal(B, E, value);
                e.addH3K27AcValue(h3k27);
            } else {
                // should never happen, actually
                System.out.printf("begin:%d, end:%d, enhancerStartPos:%d, enhancerEndPos: %d\n",
                        beginH3K27ac, endH3K27ac, enhancerStartPos, enhancerEndPos);
                if (beginH3K27ac < enhancerStartPos) {
                    System.out.println("H3K27ac-start is 5'");
                } else {
                    System.out.println("Enhancer-start is 5'");
                }
                if (endH3K27ac < enhancerEndPos) {
                    System.out.println("H3K27ac-end is 5'");
                } else {
                    System.out.println("Enhancer-end is 5'");
                }
                throw new RuntimeException("Missed case for overlap between H3K27Ac and Ehancer");
            }
        }
    }


    private void performChiSquareTest(long a, long b, long c, long d) {
        long A[][] = new long[2][2];
        A[0][0] = a;
        A[0][1] = b;
        A[1][0] = c;
        A[1][1] = d;
        double chi2 = chiSquare(A);
        double pval = chiSquareTest(A);
        System.out.printf("Zero values-CGI: %.2f%%, non CGI: %.2f%%\n", 100.0*(double)a/((double)(a+b)), 100.0 * (double)c/((double)(c+d)));
        System.out.printf("[INFO] chi2: %f, pval = %e\n", chi2, pval);
    }


    public void calculateMeanCGIvsNonCGI() {
        List<Double> cgi = new ArrayList<>();
        List<Double> noncgi = new ArrayList<>();
        int total = 0;
        int totalabovezero = 0;
        for (List<RegulatoryElement> enhs : chromosome2enhancerList.values()) {
            for (RegulatoryElement enh : enhs) {
                total++;
                double mean = enh.getMeanH3K27AcPer1000(number_of_experiments);
                if (mean > 0) {
                    totalabovezero++;
                }
                if (enh.isCpG()) {
                    cgi.add(mean);
                } else {
                    noncgi.add(mean);
                }
            }
        }
        DescriptiveStatistics cgistats = new DescriptiveStatistics();
        DescriptiveStatistics noncgistats = new DescriptiveStatistics();
        int zeroActivityCGI = 0; // count of CGI-items with ZERO h3K28ac
        int nonZeroActivityCGI = 0;
        int zeroActivityNonCGI = 0;// count of non-CGI-items with ZERO h3K28ac
        int nonZeroActivityNonCGI = 0;

        for (Double v : cgi) {
            if (v == 0.0) {
                zeroActivityCGI++;
            } else {
                cgistats.addValue(v);
                nonZeroActivityCGI++;
            }
        }
        for (Double v : noncgi){
            if (v == 0.0) {
                zeroActivityNonCGI++;
            } else {
                noncgistats.addValue(v);
                nonZeroActivityNonCGI++;
            }
        }
	int totalCGI = zeroActivityCGI + nonZeroActivityCGI;
        System.out.printf("[INFO] CGI: Zero: %d (%.1f%%), Nonzero: %d (%.1f%%)\n",
			  zeroActivityCGI,
			  (double)zeroActivityCGI*100.0/(double)(totalCGI),
			  nonZeroActivityCGI,
			  (double)nonZeroActivityCGI*100.0/(double)(totalCGI));
	int totalNonCGI = zeroActivityNonCGI + nonZeroActivityNonCGI;
        System.out.printf("[INFO] non-CGI: Zero: %d (%.1f%%), Nonzero: %d (%.1f%%)\n",
			  zeroActivityNonCGI,
			  (double)zeroActivityNonCGI*100.0/(double)(totalNonCGI),
			  nonZeroActivityNonCGI,
			  (double)nonZeroActivityNonCGI*100.0/(double)(totalNonCGI));
        performChiSquareTest(zeroActivityCGI, nonZeroActivityCGI, zeroActivityNonCGI, nonZeroActivityNonCGI);
        System.out.printf("[INFO] Analyzed %d regulatory elements (%d were above zero)\n", total, totalabovezero);
        System.out.printf("[INFO] Mean H3K27Ac (CGI): %f.\n", cgistats.getMean());
        System.out.printf("[INFO] Mean H3K27Ac (Non-CGI): %f.\n", noncgistats.getMean());
    }


    public void output_for_R(String filename) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
            for (List<RegulatoryElement> enhs : chromosome2enhancerList.values()) {
                for (RegulatoryElement enh : enhs) {
                    double mean = enh.getMeanH3K27AcPer1000(number_of_experiments);
                    if (mean==0.0) continue;
                    if (enh.isCpG()) {
                        bw.write(String.format("%s\t%f\n", "cgi", mean ));
                    } else {
                        bw.write(String.format("%s\t%f\n", "non.cgi", mean ));
                    }
                }
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }


}
