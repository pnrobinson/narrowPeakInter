# narrowPeakInter
Narrow Peak Intersection 

To run the app, we need the tss-stats-hg38e.txt and tss-stats-hg38p.txt
files. Adjust the paths as needed.
```aidl
$ mvn package
$ java -jar target/npi.jar -e tss-stats-hg38e.txt -p tss-stats-hg38p.txt
```

The R script located in the ``script`` folder can be used
to make plots with the results.