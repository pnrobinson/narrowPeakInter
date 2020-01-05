package org.jax.npi;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import org.jax.npi.analysis.ChromosomeWithEnhancers;
import org.jax.npi.data.RegulatoryElement;
import org.jax.npi.io.NarrowPeakDownloader;
import org.jax.npi.io.RegulatoryElementTssParser;
import org.jax.npi.io.TssEnhancerStatsParser;
import org.jax.npi.io.TssPromoterStatsParser;

import java.io.File;
import java.util.List;

public class NarrowPeakInter {

    @Parameter(names = {"-e","--enhancer"}, description = "path to tss-stats-hg38e.txt file", required = true)
    private String enhancerPath;

    @Parameter(names = {"-p","--promoter"}, description = "path to tss-stats-hg38p.txt file", required = true)
    private String promoterPath;







    public static void main(String [] argv) {
        NarrowPeakInter m = new NarrowPeakInter();
        try {
            JCommander.newBuilder()
                    .addObject(m)
                    .build().
                    parse(argv);
        } catch (ParameterException e) {
            e.printStackTrace();
            return;
        }
        try {
            m.run();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public NarrowPeakInter() {
    }

    private void analyzeH3K27ac(String group) {
        String tssFile;
        RegulatoryElementTssParser tssParser;
        if (group.equals("enhancer")) {
            tssFile = this.enhancerPath;
            tssParser = new TssEnhancerStatsParser(tssFile);
        } else if (group.equals("promoter")) {
            tssFile = this.promoterPath;
            tssParser = new TssPromoterStatsParser(tssFile);
        } else {
            System.err.println("[ERRROR] Did not recognize group:" + group);
            return;
        }

        List<RegulatoryElement> enhancers = tssParser.getEnhancerList();
        if (enhancers.isEmpty()) {
            throw new RuntimeException("Was not able to parse any enhancers");
        }
        ChromosomeWithEnhancers chromwe = new ChromosomeWithEnhancers(enhancers);
        File folder = new File("data");
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.getAbsolutePath().endsWith(".bed.gz")) {
                System.out.println(fileEntry.getName());
                chromwe.addDataFromBedFile(fileEntry);
            }
        }
        chromwe.calculateMeanCGIvsNonCGI();
        String outputfilename = String.format("h3k27ac-%s.txt", group);
        chromwe.output_for_R(outputfilename);
    }


    private void run() {
        NarrowPeakDownloader downloader = new NarrowPeakDownloader();
        downloader.download();
        analyzeH3K27ac("enhancer");
        System.out.println("########################################");
        analyzeH3K27ac("promoter");

    }






}
