package org.jax.npi;

import org.jax.npi.analysis.ChromosomeWithEnhancers;
import org.jax.npi.data.Enhancer;
import org.jax.npi.io.NarrowPeakDownloader;
import org.jax.npi.io.TssStatsParser;

import java.io.File;
import java.util.List;

public class Main {



    public static void main(String [] argv) {

        if (argv.length < 1) {
            System.err.println("[Usage]: java -jar npi.jar /path/to/tss-stats-hg38e.txt");
            System.exit(1);
        }
        String tssFile = argv[0];
        if (!(new File(tssFile).exists())) {
            System.err.println("[ERROR] could not find tss-stats-hg38e.txt file");
        }
        NarrowPeakDownloader downloader = new NarrowPeakDownloader();
        downloader.download();


        TssStatsParser tssParser = new TssStatsParser(tssFile);
        List<Enhancer> enhancers = tssParser.getEnhancerList();
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
        String outputfilename = "h3k27ac-enhancers.txt";
        chromwe.output_for_R(outputfilename);
    }
}
