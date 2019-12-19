package org.jax.npi.data;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class Enhancer implements Comparable<Enhancer> {

    private final String chromosome;
    private final int begin;
    private final int end;
    private final List<Double> h3k27acValues;
    private final boolean isCpG;

    public Enhancer(String chrom, int b, int e, boolean cpg) {
        this.chromosome = chrom;
        this.begin = b;
        this.end = e;
        this.h3k27acValues = new ArrayList<>();
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

    public void addH3K27AcValue(double v) {
        this.h3k27acValues.add(v);
    }

    public double getMeanH3K27Ac() {
        if (h3k27acValues.isEmpty()) return 0.0;
        return this.h3k27acValues.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
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
