package org.jax.npi.data;

import java.util.*;

public class Enhancer implements Comparable<Enhancer> {

    private final String chromosome;
    private final int begin;
    private final int end;
    private final Map<String, List<H3K27AcSignal>> experiment2H3K27map;
    private final boolean isCpG;

    public Enhancer(String chrom, int b, int e, boolean cpg) {
        this.chromosome = chrom;
        this.begin = b;
        this.end = e;
        this.experiment2H3K27map = new HashMap<>();
        this.isCpG = cpg;
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

    public void addH3K27AcValue(H3K27AcSignal signal, String experiment) {
        this.experiment2H3K27map.putIfAbsent(experiment, new ArrayList<>());
        List<H3K27AcSignal> lst = this.experiment2H3K27map.get(experiment);
        lst.add(signal);
    }

    public double getMeanH3K27AcPer1000(int expectedTotal) {
        List<Double> meanH3K27per100nt = new ArrayList<>();
        if (this.experiment2H3K27map.isEmpty())
            return 0.0;
        for (List<H3K27AcSignal> lst  : experiment2H3K27map.values()) {
            double weightedSum = 0.0;
            int totallen = 0;
            for (H3K27AcSignal h3 : lst) {
                int len = h3.getLen();
                weightedSum += len * h3.getValue();
                totallen += len;
            }
            int enhancerlen = this.end - this.begin;
            if (totallen > enhancerlen) {
                // should never happen. Sanity check
                throw new RuntimeException("BADNESS -- Total len of segments more than enhancer len");
            }
            weightedSum = weightedSum * 1000.0/enhancerlen;
            meanH3K27per100nt.add(weightedSum);
        }
        if (expectedTotal > meanH3K27per100nt.size()) {
            // in this case, some experiments are missing data for this region. Assume that
            // this means there is a zero value
            int delta = expectedTotal - meanH3K27per100nt.size();
            for (int j = 0; j < delta; j++) {
                meanH3K27per100nt.add(0.0); // adding 'fake' value for missing experiment
            }
        }
        double m = meanH3K27per100nt.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
        return m;
    }



    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Enhancer enhancer = (Enhancer) o;
        return begin == enhancer.begin &&
                end == enhancer.end &&
                Objects.equals(chromosome, enhancer.chromosome);
    }

    @Override
    public int hashCode() {
        return Objects.hash(chromosome, begin, end);
    }

    @Override
    public int compareTo(Enhancer enhancer) {
        return Integer.compare(begin, enhancer.begin);
    }
}
