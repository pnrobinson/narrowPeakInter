package org.jax.npi;

import org.jax.npi.data.Enhancer;
import org.jax.npi.io.NarrowPeakDownloader;
import org.jax.npi.io.TssStatsParser;

import java.util.List;

public class Main {



    public static void main(String [] argv) {
        NarrowPeakDownloader downloader = new NarrowPeakDownloader();
        downloader.download();
        String tssFile = "/home/peter/GIT/manuscript-enhancers2/source/tss-stats-hg38e.txt";
        TssStatsParser tssParser = new TssStatsParser(tssFile);
        List<Enhancer> enhancers = tssParser.getEnhancerList();

    }
}
