package org.jax.npi.io;

import org.jax.npi.data.RegulatoryElement;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TssEnhancerStatsParser implements RegulatoryElementTssParser {

    /** path to tss-stats-hg38e.txt file. */
    private final String pathToTssStatsFile;

    private List<RegulatoryElement> enhancerList = new ArrayList<>();

    public TssEnhancerStatsParser(String path) {
        this.pathToTssStatsFile = path;
        parse();
    }

    public List<RegulatoryElement> getEnhancerList() {
        return enhancerList;
    }


    /**
     * tss     enhancer        tags    tau.c   tau.t   dispersion      dispersion<=12  CpG.p   CpG.m
     * D       R       F       TATA    BREu    BREd    Inr     TCT     XCPE1   XCPE2   DPE     MTE     Bridge  DCE     DCE3
     * F[1] enhancer position
     * F[7] CpG.p
     * F[8] CpG.m
     * Just parse the plus strands --
     */
    private void parse() {
        int cpg_count = 0;
        int non_cpg_count = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(this.pathToTssStatsFile))) {
            String line = br.readLine(); // discard header
            while ((line = br.readLine())!= null) {
                String []F = line.split("\t");
                String tss = F[0];
                boolean isPlusStrand = tss.contains("+");
                if (! isPlusStrand) {
                    continue;
                }
                String pos = F[1];
                String cpgp = F[7];
                String cpcm = F[8];
                String []A = pos.split(":");
                if (A.length != 2) {
                    throw new RuntimeException("Bad position string: " + pos);
                }
                String chrom = A[0];
                String []B = A[1].split("-");
                if (B.length != 2) {
                    throw new RuntimeException("Bad position string: " + pos);
                }
                int start = Integer.parseInt(B[0]);
                int end = Integer.parseInt(B[1]);
                if (cpgp.equals("1") && cpcm.equals("0")) {
                    enhancerList.add(new RegulatoryElement(chrom, start, end, true));
                    cpg_count++;
                } else if (cpgp.equals("0") && cpcm.equals("1")) {
                    enhancerList.add(new RegulatoryElement(chrom, start, end, false));
                    non_cpg_count++;
                } else {
                    System.err.printf("Bad code cpcp=%s cpcm=%s\n", cpgp, cpcm);
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.printf("[INFO] Parsed %d CGI enhancers and %d non-CGI enhancers.\n", cpg_count, non_cpg_count);
    }
}
