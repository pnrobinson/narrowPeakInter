package org.jax.npi.data;


/**
 * This represents one area of H3K27ac signal in one experiment.
 */
public class H3K27AcSignal {

    public final int begin;
    public final int end;
    public final double value;

    public H3K27AcSignal(int b, int e, double v) {
        this.begin = b;
        this.end = e;
        this.value = v;
    }

    public int getBegin() {
        return begin;
    }

    public int getEnd() {
        return end;
    }

    public double getValue() {
        return value;
    }

    public int getLen() {
        // BED format no need to substract 1
        return end - begin;
    }
}
