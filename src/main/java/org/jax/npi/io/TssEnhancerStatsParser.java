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
     * (location) CpG+_anystrand  CpG+_-strand    CpG+_+strand
     * chr10:100006233-100006603       0       0       0
     * chr10:100008181-100008444       0       0       0

     */
    private void parse() {
        int cpg_count = 0;
        int non_cpg_count = 0;
        if (! this.pathToTssStatsFile.contains("cpg-hg38e.txt")) {
            throw new RuntimeException("Enhancer file must be 'cpg-hg38e.txt' and not 'tss-stats-hg38e.txt'");
        }
        try (BufferedReader br = new BufferedReader(new FileReader(this.pathToTssStatsFile))) {
            String line = br.readLine(); // discard header
            while ((line = br.readLine())!= null) {
                String []F = line.split("\t");
                if (F.length < 4) {
                    System.err.printf("[ERROR] malformed line with %d fields: %s\n", F.length, line);
                    continue;
                }
                String pos = F[0];
                String cpg = F[1];
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
                if (cpg.equals("1")) {
                    enhancerList.add(new RegulatoryElement(chrom, start, end, true));
                    cpg_count++;
                } else {
                    enhancerList.add(new RegulatoryElement(chrom, start, end, false));
                    non_cpg_count++;
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.printf("[INFO] Parsed %d CGI enhancers and %d non-CGI enhancers.\n", cpg_count, non_cpg_count);
    }
}
