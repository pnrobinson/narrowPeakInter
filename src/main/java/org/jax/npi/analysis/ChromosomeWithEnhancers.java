package org.jax.npi.analysis;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jax.npi.data.Enhancer;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class ChromosomeWithEnhancers {

    Map<String, List<Enhancer>> chromosome2enhancerList;

    public ChromosomeWithEnhancers(List<Enhancer> enhancerList) {
        chromosome2enhancerList = new HashMap<>();
        for (Enhancer e : enhancerList) {
            chromosome2enhancerList.putIfAbsent(e.getChromosome(), new ArrayList<>());
            List<Enhancer> enlst = chromosome2enhancerList.get(e.getChromosome());
            enlst.add(e);
        }
        // sort the Enhancers
        for (List<Enhancer> lst : chromosome2enhancerList.values()) {
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
    }

    /** Naive algorithm, but fast enough in practice. */
    private void addDataPoint(String chrom, int begin, int end, double value) {
        if (! this.chromosome2enhancerList.containsKey(chrom)) {
            System.err.println("Could not find chromosome " + chrom);
            return;
        }
        List<Enhancer> enhancers = this.chromosome2enhancerList.get(chrom);
        for (Enhancer e : enhancers) {
            if (begin > e.getEnd()) {
                return; // if this is the case, there is no overlap -- the enhancers are sorted
            } else if (end < e.getBegin()) {
                continue; //(begin, end) is located 5' to the enhancer
            } else if (begin >= e.getBegin() || end <= e.getEnd()){
                // if we get here, there is OVERLAP
                e.addH3K27AcValue(value);
            }
        }
    }


    public void calculateMeanCGIvsNonCGI() {
        List<Double> cgi = new ArrayList<>();
        List<Double> noncgi = new ArrayList<>();
        for (List<Enhancer> enhs : chromosome2enhancerList.values()) {
            for (Enhancer enh : enhs) {
                if (enh.isCpG()) {
                    cgi.add(enh.getMeanH3K27Ac());
                } else {
                    noncgi.add(enh.getMeanH3K27Ac());
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
        System.out.printf("[INFO] Mean H3K27Ac (CGI): %f.\n", cgistats.getMean());
        System.out.printf("[INFO] Mean H3K27Ac (Non-CGI): %f.\n", noncgistats.getMean());
    }


}
