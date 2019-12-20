package org.jax.npi.analysis;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jax.npi.data.RegulatoryElement;
import org.jax.npi.data.H3K27AcSignal;

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
    private void addDataPoint(String chrom, int begin, int end, double value) {
        if (! this.chromosome2enhancerList.containsKey(chrom)) {
            System.err.println("Could not find chromosome " + chrom);
            return;
        }
        if (begin > end) {
            throw new RuntimeException(String.format("Begin=%d and end =%d\n", begin, end));
        }
        List<RegulatoryElement> enhancers = this.chromosome2enhancerList.get(chrom);
        for (RegulatoryElement e : enhancers) {
            if (e.getBegin() > end) {
                return; // the enhancers are sorted -- the current enhancer is already beyond the end
                // of the new interval and there was no match
            } else if (end < e.getBegin()) {
                continue; //(begin, end) is located 5' to the enhancer
            } else if (end >= e.getBegin() && end <= e.getEnd()){
                // if we get here, there is OVERLAP -- end is within the enhancer
                //int len = getOverlap(begin, end, e);
                int B = begin < e.getBegin() ? e.getBegin() : begin;
                int E = end > e.getEnd() ? e.getEnd() : end;
                H3K27AcSignal h3k27 = new H3K27AcSignal(B, E, value);
                e.addH3K27AcValue(h3k27);
            } else if (begin >= e.getBegin() && begin <= e.getEnd()){
                // if we get here, there is OVERLAP -- begin (at least) is within the enhancer
                //int len = getOverlap(begin, end, e);
                int B = begin < e.getBegin() ? e.getBegin() : begin;
                int E = end > e.getEnd() ? e.getEnd() : end;
                H3K27AcSignal h3k27 = new H3K27AcSignal(B, E, value);
                e.addH3K27AcValue(h3k27);
            }
        }
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
        for (Double v : cgi) {
            cgistats.addValue(v);
        }
        for (Double v : noncgi){
            noncgistats.addValue(v);
        }
        System.out.printf("[INFO] Analyzed %d regulatory elements (%d were above zero)\n.", total, totalabovezero);
        System.out.printf("[INFO] Mean H3K27Ac (CGI): %f.\n", cgistats.getMean());
        System.out.printf("[INFO] Mean H3K27Ac (Non-CGI): %f.\n", noncgistats.getMean());
    }


    public void output_for_R(String filename) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
            for (List<RegulatoryElement> enhs : chromosome2enhancerList.values()) {
                for (RegulatoryElement enh : enhs) {
                    if (enh.isCpG()) {
                        double mean = enh.getMeanH3K27AcPer1000(number_of_experiments);
                        bw.write(String.format("%s\t%f\n", "cgi", mean ));
                    } else {
                        double mean = enh.getMeanH3K27AcPer1000(number_of_experiments);
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
