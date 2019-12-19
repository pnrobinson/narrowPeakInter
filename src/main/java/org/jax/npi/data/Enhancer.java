package org.jax.npi.data;

import java.util.Objects;

public class Enhancer implements Comparable<Enhancer> {

    private final String chromosome;
    private final int begin;
    private final int end;
    private final boolean isCpG;

    public Enhancer(String chrom, int b, int e, boolean cpg) {
        this.chromosome = chrom;
        this.begin = b;
        this.end = e;
        this.isCpG = cpg;
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
