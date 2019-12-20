package org.jax.npi.data;

import java.util.*;

/**
 * Represents a promoter or enhancer with its H3K27ac signals
 */
public class RegulatoryElement implements Comparable<RegulatoryElement> {

    private final String chromosome;
    private final int begin;
    private final int end;
    private final List<Double> means;
    private final List<Double> maxima;
    private List<H3K27AcSignal> signals;
    private final boolean isCpG;

    public RegulatoryElement(String chrom, int b, int e, boolean cpg) {
        this.chromosome = chrom;
        this.begin = b;
        this.end = e;
        if (b > e) {
            System.err.printf("[ERROR] Could not construct enhancer with begin>end.\n");
            System.err.printf("[ERROR] chrom=%s, begin=%d, end=%d CpG=%b\n", chrom, begin, end, cpg);
        }
        this.means = new ArrayList<>();
        this.maxima = new ArrayList<>();
        this.isCpG = cpg;
        signals = new ArrayList<>();
    }



    public String getChromosome() {
        return chromosome;
    }

    public int getBegin() {
        return begin;
    }

    public int getEnd() {
        return end;
    }

    public boolean isCpG() {
        return isCpG;
    }

    public void addH3K27AcValue(H3K27AcSignal signal) {
        signals.add(signal);
    }

    public void processLastExperiment() {
        // process the H3K27AcSignal from the last experiment and then reset the list
        if (this.signals.isEmpty()) {
            this.means.add(0.0);
            this.maxima.add(0.0);
            return; // no need to reset
        }
        double max = this.signals.stream().mapToDouble(H3K27AcSignal::getValue).max().orElse(1.0);
        this.maxima.add(max);
        double weightedSum = 0.0;
        int totallen = 0;
        for (H3K27AcSignal h3 : this.signals) {
            int len = h3.getLen();
            weightedSum += len * h3.getValue();
            totallen += len;
        }
        int enhancerlen = this.end - this.begin;
        if (totallen > enhancerlen) {
            System.err.println("totallen="+totallen);
            int cumulative = 0;
            for (H3K27AcSignal h3 : this.signals) {
                int len = h3.getLen();
                cumulative += len;
                System.err.printf("Len=%d (cumulative=%d).\n", len, cumulative);
            }
            // should never happen. Sanity check
            throw new RuntimeException("BADNESS -- Total len of segments more than enhancer len");
        }
        weightedSum = weightedSum * 1000.0/enhancerlen;
        means.add(weightedSum);
        this.signals = new ArrayList<>(); // reset
    }

    public double getMeanH3K27AcPer1000(int expectedTotal) {
        if (expectedTotal > means.size()) {
            // in this case, some experiments are missing data for this region. Assume that
            // this means there is a zero value
            int delta = expectedTotal - means.size();
            for (int j = 0; j < delta; j++) {
                means.add(0.0); // adding 'fake' value for missing experiment
            }
        }
        double m = means.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
        return m;
    }

    public double getMeanMaxH3K27ac(int expectedTotal) {
        if (expectedTotal > means.size()) {
            // in this case, some experiments are missing data for this region. Assume that
            // this means there is a zero value
            int delta = expectedTotal - means.size();
            for (int j = 0; j < delta; ++j) {
                maxima.add(0.0);
            }
        }
        double m = maxima.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
        return m;
    }



    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        RegulatoryElement enhancer = (RegulatoryElement) o;
        return begin == enhancer.begin &&
                end == enhancer.end &&
                Objects.equals(chromosome, enhancer.chromosome);
    }

    @Override
    public int hashCode() {
        return Objects.hash(chromosome, begin, end);
    }

    @Override
    public int compareTo(RegulatoryElement enhancer) {
        return Integer.compare(begin, enhancer.begin);
    }
}
