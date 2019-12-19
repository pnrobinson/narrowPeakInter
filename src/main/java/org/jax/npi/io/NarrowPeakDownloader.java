package org.jax.npi.io;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.css.CSSFontFaceRule;

import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * Command to download the {@code hp.obo} and {@code phenotype.hpoa} files that
 * we will need to run the LIRICAL approach.
 */
public class NarrowPeakDownloader {
    private static final Logger logger = LoggerFactory.getLogger(NarrowPeakDownloader.class);
    /** Directory to which we will download the files. */
    private final String downloadDirectory = "data";
    /**     Homo sapiens hepatocyte originated from H9 */
    String ENCFF757CYP = "https://www.encodeproject.org/files/ENCFF757CYP/@@download/ENCFF757CYP.bed.gz";
    /**     Homo sapiens neural progenitor cell originated from H9 */
    String ENCFF779WYN = "https://www.encodeproject.org/files/ENCFF779WYN/@@download/ENCFF779WYN.bed.gz";
    /**     Homo sapiens trophoblast cell originated from H1 */
    String ENCFF698NII = "https://www.encodeproject.org/files/ENCFF698NII/@@download/ENCFF698NII.bed.gz";
    /**     Homo sapiens mesendoderm originated from H1 */
    String ENCFF459UTL = "https://www.encodeproject.org/files/ENCFF459UTL/@@download/ENCFF459UTL.bed.gz";
    /**     Homo sapiens neural stem progenitor cell originated from H1 */
    String ENCFF874YBQ = "https://www.encodeproject.org/files/ENCFF874YBQ/@@download/ENCFF874YBQ.bed.gz";
    /**     Homo sapiens mesenchymal cell originated from H1 */
    String ENCFF196AMI ="https://www.encodeproject.org/files/ENCFF196AMI/@@download/ENCFF196AMI.bed.gz";
    /**     Homo sapiens endodermal cell originated from HUES64 */
    String ENCFF587KQG = "https://www.encodeproject.org/files/ENCFF587KQG/@@download/ENCFF587KQG.bed.gz";
    /**     Homo sapiens mesodermal cell originated from HUES64 */
    String ENCFF812JNL = "https://www.encodeproject.org/files/ENCFF812JNL/@@download/ENCFF812JNL.bed.gz";
    /**     Homo sapiens ectodermal cell originated from HUES64 */
    String ENCFF110UVX = "https://www.encodeproject.org/files/ENCFF110UVX/@@download/ENCFF110UVX.bed.gz";
    /**     Homo sapiens bipolar neuron originated from GM23338 treated with 0.5 μg/mL doxycycline hyclate for 4 days*/
    String ENCFF783DOC = "https://www.encodeproject.org/files/ENCFF783DOC/@@download/ENCFF783DOC.bed.gz";
    /**     Homo sapiens neuroepithelial stem cell genetically modified using stable transfection originated from H9*/
    String ENCFF168FUG ="https://www.encodeproject.org/files/ENCFF168FUG/@@download/ENCFF168FUG.bed.gz";
    /**     Homo sapiens neural cell originated from H1 */
    String ENCFF088CLP = "https://www.encodeproject.org/files/ENCFF088CLP/@@download/ENCFF088CLP.bed.gz";
    /**     Homo sapiens myotube originated from skeletal muscle myoblast */
    String ENCFF626ZXA = "https://www.encodeproject.org/files/ENCFF626ZXA/@@download/ENCFF626ZXA.bed.gz";
    public NarrowPeakDownloader(){
    }

    /**
     * Download the files unless they are already present.
     */
    public void download() {
        downloadFileIfNeeded(ENCFF757CYP);
        downloadFileIfNeeded(ENCFF779WYN);
        downloadFileIfNeeded(ENCFF698NII);
        downloadFileIfNeeded(ENCFF459UTL);
        downloadFileIfNeeded(ENCFF874YBQ);
        downloadFileIfNeeded(ENCFF196AMI);
        downloadFileIfNeeded(ENCFF587KQG);
        downloadFileIfNeeded(ENCFF812JNL);
        downloadFileIfNeeded(ENCFF110UVX);
        downloadFileIfNeeded(ENCFF783DOC);
        downloadFileIfNeeded(ENCFF168FUG);
        downloadFileIfNeeded(ENCFF088CLP);
        downloadFileIfNeeded(ENCFF626ZXA);

    }


    /**
     * From this: https://www.encodeproject.org/files/ENCFF757CYP/@@download/ENCFF757CYP.bed.gz
     * extract this data/ENCFF757CYP.bed.gz
     * @param webAddress
     * @return
     */
    private String extractLocalName(String webAddress) {
        String A[] = webAddress.split("/");
        String localname = A[A.length-1];
        return String.format("%s%s%s", downloadDirectory, File.separator, localname);
    }


    private void downloadFileIfNeeded(String webAddress) {
        String localName = extractLocalName(webAddress);
        File f = new File(localName);
        if (f.exists()) {
            System.out.println(String.format("Cowardly refusing to download %s since we found it at %s",
                    localName,
                    f.getAbsolutePath()));
            logger.trace(String.format("Cowardly refusing to download %s since we found it at %s",
                    localName,
                    f.getAbsolutePath()));
            return;
        }
        FileDownloader downloader=new FileDownloader();
        try {
            URL url = new URL(webAddress);
            logger.debug("Created url from "+webAddress+": "+url.toString());
            downloader.copyURLToFile(url, new File(f.getAbsolutePath()));
        } catch (MalformedURLException e) {
            logger.error(String.format("Malformed URL for %s [%s]",localName, webAddress));
            logger.error(e.getMessage());
        } catch (FileDownloadException e) {
            logger.error(String.format("Error downloading %s from %s" ,localName, webAddress));
            logger.error(e.getMessage());
        }
        System.out.println("[INFO] Downloaded " + localName);
    }





}
