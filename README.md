# 850K-EPIC-analysis-pipeline-FRANSQUET

DISCLAIMER: This is to be used as a guide only! 
Much of my code may not be pretty, some methods may be outdated, some packages may need to be updated or not currently work, and some analysis methods or statistics used may be incorrect or not appropriate for your dataset. Please feel free to reach out if there is something terribly wrong. 

This is a pipeline I have put together based off the amazing work by Jovana Maksimovic, Belinda Phipson and Alicia Oshlack
https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
Work was carried out using EPIC EWAS during my PhD (2017-2020). I used 850k data from blood methylation but much of the methods could also be used for other tissues or 450k datasets.


Methods within the pipeline include:
- Reading in the data
- Quality control
  -Normalisation
  -Some basic data visualization
  -Checking sample sex
  -Removing problematic probes
- Creating M and Beta data sets
- Cell estimation
- Principal Component Analysis (PCA)
- Epigenetic Age estimation using Horvath’s Online Calculator
- Analysis using Limma/'cate' for categorical or continuous regression models
- Candidate gene analysis
- Differentially Methylated Region (DMR) analysis using 'DMRcate'
- Gene pathway analysis using 'gometh'
