#################################################### Version 4.1 - 14/NOV/2022 ####################################################
## Analysis Pipeline for Illumina EPIC (850K) Methylation Array
## This is based off 450k pipeline by Jovana Maksimovic*, Belinda Phipson and Alicia Oshlack
## https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
## Have a read as there are some things that are done there that are not in this pipeline

## Also this code uses updated annotation "IlluminaHumanMethylationEPICanno.ilm10b4.hg19" which caused some issues with the pipeline because
## there were different amount of probes between what is on the chip and what is in the annotation (mset had a different amount of probes compared to rgset at one stage)
## I wrote a probe list work which fixes it, but can still use "IlluminaHumanMethylationEPICanno.ilm10b2.hg19" if you dont want to use the work around.


################################################  Some simple handy Pete notes   #################################################

## Some things to keep in mind: 

## Some things take ages to run, and depend on the speed of the computer you use. Much of this code, depending on how many samples, will mean you cannot run on even a great home PC
## and will need access to an RStudio server which are linked up to super computers. Some of these servers are used simultaneously by others which will slow your work down.

## Every data set is different and will have different needs and issues, so you will likely never just be able to cut and paste and have it work
## This pipeline needs to be interpreted/translated for each individual dataset/server.

## There are lots of things that can go wrong in R and with programming in general, sometimes the errors don't make sense, some errors come up but it doesnt mean anything for the progression of the pipeline.
## Sometimes packages don't load from particular sources, R needs updating, file types aren't as they should be etc. There are too many to list here.
## That said, I spend a lot of time googling things. There is usually a fix for everything, just need to ask the right question and translate and interpret what you find.

## if you have never used R, you can use ctrl+ENTER to run a line of code, or highlight several lines then use ctrl+ENTER to run the selection.
## you can find information on a package or function by using two questionmarks in front of its name.
## eg ## ??plotMDS


## If you want a fresh start you can Wipe the Global Environment using ## rm(list=ls()) ##
## ctrl+l will clear the console
## ctrl+shift+F10 will restart R, dumping everything in RAM

## you can save a data frame (named DF for examples) using for example ## save(DF,file="c:/mydir/myDF.Rdata") ##
## and then remove it from your environment by using ## rm(DF) ##

## To load again use ## load("c:/mydir/myDF.Rdata") ## or by finding it in your folder and clicking on it to load it back in
## This will load a DF rather than having to build it again

## Keep in mind your file location might be different depending on how your RStudio is set up, could be on a hard drive or on a server.
## so could be something more like ## load("/mydir/myDF.Rdata") ##

## Sometimes data is in a matrix but needs to be a data frame and vice versa, but the errors are often and ambiguous
## to check the R file type and the types of data within you can use ## str(dataframe_name) ##

## To remove a column from a dataframe (using "Varname" as a column), simply write # DF$Varname = NULL ##

## To get rid of a row based on a condition, for example, get rid of a row based on the variable name 256
## DF<-DF[!(DF$Varname=="256"),]

## To write a data frame into a "*.csv" file
## write.csv(DF,'DF.csv')


####################################################  Editing the Raw EWAS data and SampleSheet ##################################################### 
## Extract all compressed (zipped) files into one singular folder so all the unzipped data is all in one place. I usually put this in a folder named 'EWAS' in whatever RStudio server i am using.

## Now you need to do is fix up your "SampleSheet" that you get from whoever has run the assay (AGRF in my case), and turn into a 'targets' file, saved as a csv.
## I do all the edits manually using excel:

## 1. Make a copy of the sample sheet, and be sure you do not edit the original, keep the original somewhere safe.
## 2. Open the copy of the sample sheet, remove all the data not needed which is everything above the main data (headings including Sample_Name, Sample_Well, Sentrix_ID etc.)
## 3. VERY IMPORTANT: Make sure to specify Sentrix_ID column as 'number' format with no decimals, or excel will make it 'scientific' format losing all the data which tells R what sample is where
## 4. Add any participant data to the file that you will need, such as age, sex, outcome and possible confounding variables.
## 5. Once completed, save as a "CSV (Comma Delimited) (*.csv)" file named 'targets.csv', as R will want to read it in this format. Do not use any other CSV format (for example MACINTOSH) or it will not work.
## 6. Copy this 'targets.csv' into the same folder with all the extracted raw EWAS data.

## You might find you want to add or change things, along the way you can add using R and save for export to use again in excel or change in excel and start over.
## For example i use a grouped ID where i have the sample ID + case or control e.g. 24_C for a control and 25_D for a case of dementia.
## Some notes: If spacing is needed in either the table data or headings use "_", a period "." makes R think its a weird number, and "-" changes things to a date.
## R will not like 0 in some cases when using binary coding (for example male may be coded as 0 and female coded as 1), so I changed mine (0's to 1's) and (1's to 2's)



############################################## Installing and Loading needed packages #############################################
## There are different methods for loading in packages, most of the time you need to be online, if using a secure interface like "SafeHaven" at monash you might need to get 
## Someone from Safehaven to install for you.
## Here you might get errors (e.g. Error: Bioconductor version '3.12' requires R version '4.0'; see https://bioconductor.org/install) 
## or need updates, but generally you can use ### BiocManager::install("your_package_here") ###

## Can check the version of BiocManager by using:
BiocManager::version()
## you may need to update BiocManager depending on R version, can find info here ## https://bioconductor.org/install
## Also may need to use ## setRepositories() ## to be able to access other sources to download packages

## Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
## as of writing this version 3.14 is the newest version for R 4.1.1


## Sometimes a package you need isn't available for the newer version so you might have to mix and match which gets tricky.
## You can check package versions using ## packageVersion('your_package_here') ##
## Here is the list of all the packages I have used at some stage or another, you may not need them all, and sometimes they can take a while to install.
## I suggest installing one package at a time, as many of these packages need bits and pieces from other packages and will often fail to install if doing it all in one shot.
## Some packages may install while loading other packages too.

BiocManager::install("limma")
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19") ## Might also need ## BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19" ## depending on process.
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("RColorBrewer")
BiocManager::install("lumi")
BiocManager::install("missMethyl")
BiocManager::install("matrixStats")
BiocManager::install("Gviz")
BiocManager::install("DMRcate")
BiocManager::install("stringr")
BiocManager::install("WGCNA")
BiocManager::install("AnnotationDbi")
BiocManager::install("impute")
BiocManager::install("GO.db")
BiocManager::install("preprocessCore")
BiocManager::install("ruv")
BiocManager::install("FlowSorted.Blood.EPIC")
BiocManager::install("ExperimentHub")
BiocManager::install("ggplot2")
BiocManager::install("cate")
BiocManager::install("sva")
BiocManager::install("plyr")
BiocManager::install("dplyr")
BiocManager::install("forcats")
aBiocManager::install("ggpubr")


## Now need to load each library, i also do this one at a time because of crossover
## Also if you have to many big libraries loaded you can overload the amount of RAM, talk to system admin (Jerico for me) if you need more
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(lumi)
library(missMethyl)
library(matrixStats)
library(Gviz)
library(DMRcate)
library(stringr)
library(WGCNA)
library(AnnotationDbi)
library(impute)
library(GO.db)
library(preprocessCore)
library(ruv)
library(FlowSorted.Blood.EPIC)
library(ExperimentHub)
library(ggplot2)
library(cate)
library(sva)
library(plyr)
library(dplyr)
library(forcats)
library(ggpubr)




############################################## Loading in sample and methylation files ##############################################

## Set working directory, and set base directory to where the EWAS data and 'targets.csv' file is on the Hard Drive/server,
## Can check using 'list.files(baseDir)' to look at the files in the directory and make sure you are in the right spot
setwd ("/misc/mcclr-lrcbiodata/ATP_EWAS/Raw Extracted All")
baseDir="ewas"
list.files(baseDir)

## Next load in the Illumina Epic Probe annotation into a file named 'ann.EPIC'
## It's basically all the information on each probe, including name, location, SNP status and gene name etc.
ann.EPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
## You can use the ## head() ## code to see what the annodation looks like at a glance.
head(ann.EPIC)

## You need to premake folders to save into such as the ones i have below (/processing/R Data/Core/)
## I also make folders for QC, CSV/excel files, and graphics outputs etc.
## In this code i have certain files in certain places, change based on your own preferences.

## Specifically I make a 'Core' folder where I keep main files i will use, so i can:
## 1) bring them back without needing to rebuild them and 2) i can remove them from the RStudio environment when not in use which will save on the system's RAM.
## Most of the files in this first part will be saved as core files, and take up gigs and gigs worth of space
## I usually put the date in the name, sometimes rebuilds are needed and this is a way to serialize them. I also put in 'B4' as this is the annotation i used from illumina.
## Save the ann.EPIC you just made as an '*.Rdata' file so you can reload straight back in when needed.
save(ann.EPIC,file="~/processing/R Data/Core/ann.EPICB4.04OCT2021.Rdata")

## Now need to load in that targets.csv file made previously from the SampleSheet.xls file in the base directory##
targets = read.metharray.sheet(baseDir, "targets.csv")
## Can use the View() code to look at your targets file, make sure its all looking ok
View(targets)
## and str() code to see what type of data is in the data frame, sometimes importing might put what you thought is a number as a char so this is a good way to see what was done.
str(targets)

## Now need to make a data set of the methylation intensities from the raw data.
## This can take ages depending on the computer you are using, so be patient
## From the Rstudio server though Monash, it is currently around 20 samples per 1 min or so, but some systems it has gone as slow as 6 samples a minute.
rgset = read.metharray.exp(targets=targets, force=TRUE)
head(rgset)

## You want to now name the samples in the rgset, the same name as what you will be using along the way.
## In my case I made a "group" ID as GID which i will be using as the identifier.
sampleNames(rgset) = targets$GID

## Now save the files you just made, also making a new targets.csv file you can work with in the future if needed.
write.csv(targets,file="~/processing/R Data/CSV Files/targets04OCT2021.csv")
save(targets,file="~/processing/R Data/Core/targets04OCT2021.Rdata")
save(rgset,file="~/processing/R Data/Core/rgset04OCT2021.Rdata")


############################################## Basic Quality control ##############################################
## Next you create a density plot to see the distribution of the data
## If you have ran more than 8 chips (so 64 samples), and you want to break down by chip for colours, you need to change the colour scheme or it will only show 8 chips.

## load what ever colours you want into the R Colour palette, this is just a bunch of random ones i chose, this adds 47 colours so it will be ok as long as there is not more than 47+8 chips.
palette("default")
cc <- palette()
palette(c(cc,"red","orange", "wheat4","aquamarine3", "violet", 
          "cyan4", "darkblue", "brown3", "darkorchid4", "goldenrod1", "darkturquoise", "deeppink4", "green", "green4",
          "lightcoral", "hotpink", "lightslateblue", "navy",
          "orange3", "orangered2", "red4", "royalblue", "palevioletred",
          "sandybrown", "plum", "seagreen3", "peru", "springgreen4", "slategray4",
          "slateblue2", "tomato3", "violetred", "wheat3", "sienna", "tan3",
          "pink3", "peachpuff2", "salmon", "paleturquoise3"))
cc <- palette()

##see what colours you have on a plot
plot(1:51, 1:51, col=1:51, pch=19, cex=3, xlab="", ylab="")

## Now to Create a density plot (bi-modal methylation) of raw rgset data, colour seperated by chip number (which is named 'Slide' in the original sample sheet)
## Can also colour by different variables depending on how you want to visualize data. 
## Just change the identifier in the "sampGroups = " term.
## For example you could have 'sampGroups = targets$sex', which would colour by male and female.
pdf(file="~/processing/R Data/QC/RGSET Density plot 04OCT2021 (Coloured by Chip number).pdf", height=8, width=12)
densityPlot(rgset, sampGroups = targets$Slide, pal = palette(), legend = FALSE)
dev.off()
## The output is the raw distribuion of methylation between 0 and 1. Sometimes the raw looks great and sometimes it doesn't, which is why we normalise down the track.

## Now we want to use detectionP to helps find failed samples/probes.
## It calculates detection P-values, this function identifies failed positions defined as
## "both the methylated and unmethylated channel reporting background signal levels", "p values over 0.01 should not be trusted"
## more info here https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/detectionP

## i use 'options(scipen=999)' because if you don't do this there will be scientific notation in the tables! can change back to smaller number if needed.
## you may have to change this options(scipen=999) depending on what you are doing and how much scientific notation you want in your datasets.
options(scipen=999)
## This can also take a while
detp = detectionP(rgset)
head(detp)
save(detp,file="~/processing/R Data/Core/detp04OCT2021.Rdata")

## Make a bar graph of each samples detp value, and visualise where you want to set a cut off depending on data, minimum cut of as discussed by the creators is >0.01.
## Note that this code is specific for my data set. you can change the colour group using col=, the y axis limit with ylim=, and where to put the cutoff line you want visualied
## within the abline(h=0.xx). Keep this cutoff number in mind for publication if needed.
pdf(file="~/processing/R Data/QC/Detection P plot all final samples04OCT2021.pdf", height=8, width=12)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,1))
barplot(colMeans(detp), col=cc[factor(targets$Pregnancy_PTSD)], las=2, ylim=c(0,0.003),
        cex.names=0.4,ylab="Mean detection p-values")
abline(h=0.0015,col="red")
legend("topleft", legend=levels(factor(targets$Pregnancy_PTSD)), fill=cc,
       bg="white")
dev.off()

## Can also make a QC report, which can help weed out some of the samples which should be removed because they might have poor data.
## Here you can colour again how you want using 'sampGroups=', might be handy to choose chip to see if there are any noticeably bad batches, but here i have done by PTSD status.
## I have never really used this because it takes a bit to understand, and there has been easier ways to pick out dud samples.
## This can take a few minutes too.
qcReport(rgset,sampNames=targets$GID, sampGroups=targets$Pregnancy_PTSD,pdf="~/processing/R Data/QC/QC report by Pregnancy_PTSD 04OCT2021.pdf")
dev.off()



############################################## Normalisation ##############################################
## Now we need to creates a 'MethylSet' object from the raw rgset data, converting the red/green channel for an Illumina methylation array into methylation signal, first without using any normalisation.
mset.raw = preprocessRaw(rgset)
save(mset.raw,file="~/processing/R Data/Core/mset.raw.04OCT2021.Rdata")

## You can use MDS plots to examine major sources of variation in the data, (mdsPlot is a minfi function, which requires an MethylSet which we just made)
## You can group (colour), and name samples depending on how you want to explore the variation, for example here i want to see if colours separate by the exposure (pregnancy PTSD), and also label by gender
## In this particular plot it is fully separated by sex, with male on one side and female on the other.
## If you notice that there is some sex discrepancy you can put the ID in 'sampNames=' so you can see which sample it is that is out.
pdf(file="~/processing/R Data/QC/MDS plot, msetraw, grouped by Pregnancy_PTSD, Gender_Child label.pdf", height=8, width=12)
mdsPlot(mset.raw, sampGroups = targets$Pregnancy_PTSD, sampNames=targets$Gender_Child, pal = palette(), legendPos= NULL)
dev.off()
## There may be other sources of variation based on your data set, for example, in one of mine because there were several runs at AGRF
## They all seperated and it could easily be seen in the mdsPlot.

## There are many normalisation methods, as far as i have seen subset quantile normalisation (SQN) and Subset-quantile Within Array Normalization (SWAN) are the most commonly used. Both have specific uses, but both work well.
## Some papers have shown that method of normalisation have no significant effect on analysis outcomes.
## I use a few methods and compare to see which "normalises" better, looking for something more uniform.
## Read SWAN: Subset-quantile Within Array Normalization for Illumina Infinium HumanMethylation450 BeadChips - Maksimovic et al. Genome Biology 2012

## For SQN (this took me about 2 min for 120 samples)
mset.SQN <- preprocessQuantile(rgset)

## For SWAN (this also took me about 2 min for 120 samples)
mset.SWAN = preprocessSWAN(rgSet = rgset, mSet = mset.raw, verbose=TRUE)

## Save the normalised data sets
save(mset.SWAN,file="~/processing/R Data/Core/mset.SWAN.04OCT2021.Rdata")
save(mset.SQN,file="~/processing/R Data/Core/mset.SQN.04OCT2021.Rdata")

## Make another density plot (bi-modal plot) to compare normalised data to raw data to see which normalises best
## When using SQN normalisation it changes the data from a methylset to a GenomicRatioSet, so we need to plot by telling it to create the beta values
pdf(file="~/processing/R Data/QC/Density plot, raw vs swan vs sqn, all samples 04OCT2021.pdf", height=6, width=18)
par(mfrow=c(1,3))
densityPlot(mset.raw, sampGroups = targets$Pregnancy_PTSD,main="RAW", pal = palette(), legend = FALSE)
densityPlot(mset.SWAN, sampGroups = targets$Pregnancy_PTSD,main="SWAN", pal = palette(), legend = FALSE)
densityPlot(getBeta(mset.SQN), sampGroups = targets$Pregnancy_PTSD,main="SQN", pal = palette(), legend = FALSE)
dev.off()


############################################## Visualizing Data ###########################################
## We can visualize the data similarly as we did before and compare each normalisation method and compare principle components (weighted variation of data) against each other. 
## This method uses limma function plotMDS,which makes a PCA or PCoA plot. See function details for more information ## ??plotMDS ##
## It may look quite similar as the main source of variation was sex and sex probes have not been removed yet
## Here we compare principle components 1 and 2, identifying data by PTSD and sex.

pdf(file="~/processing/R Data/QC/MDSPlot, all msets pre sex removal, Pregnancy_PTSD and Gender_Child 04OCT2021.pdf", height=15, width=10)
par(mfrow=c(3,2))
plotMDS(getM(mset.raw), top=1000, gene.selection="common", main="RAW by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)])
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.raw), top=1000, gene.selection="common", main="RAW by Gender_Child", col=cc[factor(targets$Gender_Child)])
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SQN), top=1000, gene.selection="common", main="SQN by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)])
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SQN), top=1000, gene.selection="common", main="SQN by Gender_Child", col=cc[factor(targets$Gender_Child)])
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SWAN), top=1000, gene.selection="common", main="SWAN by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)])
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SWAN), top=1000, gene.selection="common", main="SWAN by Gender_Child", col=cc[factor(targets$Gender_Child)])
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
dev.off()

##Can also compare different principle components as below (1-4), still identified by PTSD and sex.
## I have made seperate files for RAW, SQN and SWAN data
## Higher PC dimentions Pregnancy_PTSD RAW
pdf(file="~/processing/R Data/QC/MDS, Higher PC dimentions Pregnancy_PTSD and Gender_Child RAW 04OCT2021.pdf", height=10, width=15)
par(mfrow=c(2,3))
plotMDS(getM(mset.raw), top=1000, gene.selection="common", main="RAW PC1 vs PC3 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.raw), top=1000, gene.selection="common", main="RAW PC2 vs PC3 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(2,3))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.raw), top=1000, gene.selection="common", main="RAW PC3 vs PC4 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.raw), top=1000, gene.selection="common", main="RAW PC1 vs PC3 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.raw), top=1000, gene.selection="common", main="RAW PC2 vs PC3 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(2,3))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.raw), top=1000, gene.selection="common", main="RAW PC3 vs PC4 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
dev.off()

## Higher PC dimentions Pregnancy_PTSD SQN
pdf(file="~/processing/R Data/QC/MDS, Higher PC dimentions Pregnancy_PTSD and Gender_Child SQN 04OCT2021.pdf", height=10, width=15)
par(mfrow=c(2,3))
plotMDS(getM(mset.SQN), top=1000, gene.selection="common", main="SQN PC1 vs PC3 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SQN), top=1000, gene.selection="common", main="SQN PC2 vs PC3 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(2,3))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SQN), top=1000, gene.selection="common", main="SQN PC3 vs PC4 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SQN), top=1000, gene.selection="common", main="SQN PC1 vs PC3 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SQN), top=1000, gene.selection="common", main="SQN PC2 vs PC3 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(2,3))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SQN), top=1000, gene.selection="common", main="SQN PC3 vs PC4 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
dev.off()

## Higher PC dimentions Pregnancy_PTSD SWAN
pdf(file="~/processing/R Data/QC/MDS, Higher PC dimentions Pregnancy_PTSD and Gender_Child SWAN 04OCT2021.pdf", height=10, width=15)
par(mfrow=c(2,3))
plotMDS(getM(mset.SWAN), top=1000, gene.selection="common", main="SWAN PC1 vs PC3 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SWAN), top=1000, gene.selection="common", main="SWAN PC2 vs PC3 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(2,3))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SWAN), top=1000, gene.selection="common", main="SWAN PC3 vs PC4 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SWAN), top=1000, gene.selection="common", main="SWAN PC1 vs PC3 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SWAN), top=1000, gene.selection="common", main="SWAN PC2 vs PC3 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(2,3))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(mset.SWAN), top=1000, gene.selection="common", main="SWAN PC3 vs PC4 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
dev.off()


## After looking at the data and at which density plot looks more uniform and make a decision which nomalisation method do use
## i am going with SQN this time, mentioned in "Maksimovic 2017 bioconductor pipeline for 450K" and "Maksimovic RUV 2015" from Jane Loke
gmset = mapToGenome(mset.SQN)
save(gmset,file="~/processing/R Data/Core/gmset.SQN.04OCT2021.Rdata")
## if you look at the GM set in the Environment window it should tell you how many probes or 'elements' in the data set
## note this down (in my case it is 865859 probes) so you can report it in manuscripts as to how many probes you started with and how many you used in analsis after removing probes below
## Sometimes it bugs out and just says something like "object with Null pointer". I dont know what it means but it still works in the pipeline.
## if this is the case you can use the following code which tells you the number of rows
nrow (gmset) ##865859
## You can use this every time you remove probes to see how many you are removing



########################################### Sample Sex Check ##############################################
## You can plot sex to make sure it all groups together
pdf(file="~/processing/R Data/QC/SexCheck 04OCT2021.pdf", height=10, width=15)
plotSex(gmset, id = ifelse(targets$Gender_Child == 1, 'M', 'F'))
dev.off()

## I have also written this code, which can list off which samples are discordant with the reported sex.
sex.pred=getSex(gmset)
## In my targets, GID is col 2 and Gender_Child is col 19
sextable<- targets[,c(2,19)]
sextable$Sex <- getSex(mset.SQN)$predictedSex
sextable$ok <- ifelse(sextable$Gender_Child=="0"&sextable$Sex=="M", "OK_M",
                      ifelse(sextable$Gender_Child=="1"&sextable$Sex=="F", "OK_F", 
                             ifelse(sextable$Gender_Child=="0"&sextable$Sex=="F", "FAIL",
                                    ifelse(sextable$Gender_Child=="1"&sextable$Sex=="M","FAIL", NA))))

View(sextable)

## Now check if any and how many samples failed to match predicted sex with recorded sex (Gender_Child)
length(which(sextable$ok =="FAIL"))
## if you get "[1] 0" or intger(0)" it means none have failed!

## if you get a number then you need to get a list of the samples and figure out why
## This will create a table then a list of have failed to be concordant between predicted and recorded Gender_Child
sexdiscordant <- subset(sextable, ok=="FAIL")
sexdiscordantlist <- sexdiscordant$ID


########################################### Cleaning Data for Analysis #########################################

###########################################   Removing Failed Probes!   #########################################
## Remove any probes that have failed in one or more samples and start building a clean probe list with probes removed
## I make a data set 'gmset.pr' (probes removed)
## Remember in the publication for detp "p values over 0.01 should not be trusted" so that will be the cutoff
keep = rowSums(detp< 0.01) == ncol(gmset)
## gmset.pr = gmset[keep==TRUE,]
## If the line of code above doesn't work, and you get the error "Error: subscript is a logical vector with out-of-bounds TRUE values"
## It means that there is a different amount of probes between the gmset and what is listed in detp (represented below in mset.SQN.
## This still hasn't been fixed as of OCTOBER 2021, so i have written a work around as follows:

## First you have to create a function that creates lists of group duplicates, highlight and run this whole thing from the first line to the '}':

dupsBetweenGroups <- function (df, idcol) {
  #df: the data frame
  #idcol: the column which identifies the group each row belongs to
  
  # Get the data columns to use for finding matches
  datacols <- setdiff(names(df), idcol)
  
  # Sort by idcol, then datacols. Save order so we can undo the sorting later.
  sortorder <- do.call(order, df)
  df <- df[sortorder,]
  
  # Find duplicates within each id group (first copy not marked)
  dupWithin <- duplicated(df)
  
  # With duplicates within each group filtered out, find duplicates between groups. 
  # Need to scan up and down with duplicated() because first copy is not marked.
  dupBetween = rep(NA, nrow(df))
  dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols])
  dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols], fromLast=TRUE) | dupBetween[!dupWithin]
  
  # ============= Replace NA's with previous non-NA value ==============
  # This is why we sorted earlier - it was necessary to do this part efficiently
  
  # Get indexes of non-NA's
  goodIdx <- !is.na(dupBetween)
  
  # These are the non-NA values from x only
  # Add a leading NA for later use when we index into this vector
  goodVals <- c(NA, dupBetween[goodIdx])
  
  # Fill the indices of the output vector with the indices pulled from
  # these offsets of goodVals. Add 1 to avoid indexing to zero.
  fillIdx <- cumsum(goodIdx)+1
  
  # The original vector, now with gaps filled
  dupBetween <- goodVals[fillIdx]
  
  # Undo the original sort
  dupBetween[sortorder] <- dupBetween
  
  # Return the vector of which entries are duplicated across groups
  return(dupBetween)
}

## Now you have to generate and compare lists and remove the probes that is causing the issue.
## This might be a bit rudamentary as i wrote it in 2018 but it still works
## To get probelist from gmset
gmsetprobes <-as.data.frame(gmset@rowRanges@ranges@NAMES)
colnames(gmsetprobes)[which(names(gmsetprobes) == "gmset@rowRanges@ranges@NAMES")] <- "probes"
rownames(gmsetprobes) <- gmsetprobes$gmsetprobes
gmsetprobes$code <- "gm"
View(gmsetprobes)

## to get mset probes
msetprobes <- as.data.frame (mset.raw@NAMES)
colnames(msetprobes)[which(names(msetprobes) == "mset.raw@NAMES")] <- "probes"
rownames(msetprobes) <- msetprobes$msetprobes
msetprobes$code <- "m"
View(msetprobes)

#Stick them both together and make a duplicate identifying column
probelist <- rbind(gmsetprobes,msetprobes)                    
dupRows <- dupsBetweenGroups(probelist, "code")
probedup <-cbind(probelist, dup=dupRows)
probedupuniq <- cbind(probedup, unique=!dupRows)

uniqeprobes <- subset(probedupuniq, unique==TRUE)
rownames(uniqeprobes) <- c()
keep.probes <- gmsetprobes [,1, drop=FALSE]
dropedprobelist <- uniqeprobes[,1, drop=FALSE]

## You can save the lists for reference.
save(dropedprobelist,file="~/processing/R Data/QC/droppedlegacyprobes.04OCT2021.Rdata")
write.csv(dropedprobelist,file="~/processing/R Data/CSV/Dropped probes from mapToGenome.csv")

##Now we can remove the probes from the detp data
remove = dropedprobelist$probes
detpb4 = detp[!row.names(detp)%in%remove,]

## And pick up where we were before, removing the failed probes and note it down
keepb4 = rowSums(detpb4< 0.01) == ncol(gmset) 
gmset.pr = gmset[keepb4==TRUE,]
nrow (gmset.pr) ##847977, so this removed 17,8825 probes including failed probes and probes that werent in the gmset data compared to the mset data




###########################################   Removing SNP Probes   #########################################

##Remove probes with SNPs at CpG or SBE site with specified MAF (minor allele frequency, i chose a strict 0), you can also remove "Probe" snps but i havent found an example of this in any of the literature
## Also they do not remove "Probe" snps in the original pipeline that this is somewhat adapted from:
## https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
gmset.pr <- addSnpInfo(gmset.pr)
gmset.pr <- dropLociWithSnps(gmset.pr, snps=c("SBE","CpG"), maf=0)
nrow (gmset.pr) ##818205, so this removed 29,772 probes which are near or on SnPs




###########################################   Removing sex Probes   #########################################

#If your data includes males and females, remove the sex chromosomes
autosomes = !(featureNames(gmset.pr) %in% ann.EPIC$Name[ann.EPIC$chr %in% c("chrX","chrY")])
gmset.pr = gmset.pr[autosomes,]
nrow (gmset.pr) ##799376, so this removed 18,829 probes which are on the sex chromosomes




###########################################   Removing Cross Reactive Probes   #########################################
## Exclude cross reactive probes
## I got this list from from the first Supplementary Materials file from "Pidsley et al 2016 - Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling"
## It includes 43253 cross reactive probes
## First you need to get the file, put the heading "TargetID" above the probes and save it as 'xreactEPIC_pidsley.csv', then put it in the EWAS folder with all the other raw files
## then you can run this
cross.react = read.csv("~/ewas/xreactEPIC_pidsley.csv", stringsAsFactors=FALSE)
no.xreact = !(featureNames(gmset.pr) %in% cross.react$TargetID) 
gmset.pr = gmset.pr[no.xreact,]
## When you run this make sure to keep an eye out as there may be a more updated list.
nrow (gmset.pr) ##760590, so this removed 38,786 probes which were cross reactive (that hadn't already been removed from above).


## Now the data has been cleaned with all problematic probes removed. 
## Make sure you save this file so it is easy to come back to
save(gmset.pr,file="~/processing/R Data/Core/gmset.pr.04OCT2021.probes.removed.Rdata")

## You can rerun the PCA comparisons from above using this new data set to see if there is still a major seperation of sex or any other data
## With any luck you will still get a seperation by whatever exposure/outcome you are looking at, but most of the time it is just random (which is also ok).
pdf(file="~/processing/R Data/QC/MDS plots, SQN normalisation_ all probes removed, by Pregnancy_PTSD and Gender_Child 04OCT2021.pdf", height=10, width=20)
par(mfrow=c(2,4))
plotMDS(getM(gmset.pr), top=1000, gene.selection="common", main="Final GMset PC1 vs PC2 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(1,2))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(gmset.pr), top=1000, gene.selection="common", main="Final GMset PC1 vs PC3 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(gmset.pr), top=1000, gene.selection="common", main="Final GMset PC2 vs PC3 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(2,3))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(gmset.pr), top=1000, gene.selection="common", main="Final GMset PC3 vs PC4 by Pregnancy_PTSD", col=cc[factor(targets$Pregnancy_PTSD)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Pregnancy_PTSD)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(gmset.pr), top=1000, gene.selection="common", main="Final GMset PC1 vs PC2 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(1,2))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(gmset.pr), top=1000, gene.selection="common", main="Final GMset PC1 vs PC3 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(gmset.pr), top=1000, gene.selection="common", main="Final GMset PC2 vs PC3 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(2,3))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
plotMDS(getM(gmset.pr), top=1000, gene.selection="common", main="Final GMset PC3 vs PC4 by Gender_Child", col=cc[factor(targets$Gender_Child)], dim=c(3,4))
legend("top", legend=levels(factor(targets$Gender_Child)), text.col=cc, bg="white", cex=0.7)
dev.off()


## You might find that at this stage that some particular samples are no good or not right, you can remove them from the original targets set and start over so they are not normalized together (which is what i do)
## Or you can remove them manually using R, which i avoid just incase it is doing something i cannot see (plus i dont want bad data to pull the normalisation around)

## Also, from here, you can remove some of the other data you wont use to free up some space (you can do this at any time if you know you wont use it).
## Make sure you have saved what you need before removing things!!

## You can either delete all using:
rm(list=ls()) 
##and load back in the 'ann.EPIC', 'gmset.pr', rgset (for cell type estimation if you have blood data) and 'targets' files

## or remove everything tediously one by one using rm()
rm(mset.raw)
rm(mset.SWAN) ## etc



############################################ Making M and Beta data for Analysis and Interpretation #################################################

## Calculate M values and makes new data frame, note that dataframes are case sensitive, so if you use 'M' instead of 'm' in later analysis it will not work.
## M values are used for analysis and are calculated by = log (Methylated/Unmethylated)
m <- getM(gmset.pr)
head(m[,1:5]) ## just to check there is data there
save(m,file="~/processing/R Data/Core/m.04OCT2021.Rdata")

## Make a Beta table, beta values are used for interpretation as it is a value between 0 and 1 methylation (with 0 being no methylation and 1 being 100% methylation)
## this is calculated by = Methylated / Methylated+Unmethylated+100
beta <- getBeta(gmset.pr)
head(beta[,1:5]) ## just to check there is data there
save(beta,file="~/processing/R Data/Core/beta.04OCT2021.Rdata")

## Merge probe annotation data frame with beta values, and call it 'ann.beta'
ann.beta <- merge(beta, ann.EPIC, by = "row.names")
save(ann.beta,file="~/processing/R Data/Core/ann.beta.04OCT2021.Rdata")

## If you have a categorical case control design and want to find the difference between both groups at particular probes you can do the following
## I am using a 'status' variable that i made to make things easy for myself:
Control=beta[,targets$status == "C"]
Case=beta[,targets$status == "P"]
deltabeta=rowMeans(Case) - rowMeans(Control)
dbeta=data.frame(ID=names(deltabeta), deltabeta=deltabeta)
save(deltabeta,file="~/processing/R Data/Core/deltabeta.04OCT2021.Rdata")
## Keep this as further on you can add it to your analysis output.


## You can also visualise the M data in a cluster plot, to see if there are any further obvious out liars (or interesting similarities)
## This plot took ~2 min on the Monash RStudio server for 120 samples.
pdf(file="~/processing/R Data/QC/Clusterplot, M data, 04OCT2021.pdf", height=8, width=12)
par(cex=0.5)
plotSampleRelation(m, cv.Th = 0.01, method=c("cluster"))
dev.off()


######################################################## CELL ESTIMATION #########################################################
## I use the method from https://bioconductor.org/packages/release/data/experiment/html/FlowSorted.Blood.EPIC.html
## You need the internet for doing this section as it uses online databases, so cannot be done in places like Monash SafeHaven
## BiocManager::install("FlowSorted.Blood.EPIC", version = "devel")
## Also as of version 1.17.2 there is a change to the default location
## https://bioconductor.org/packages/devel/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html#default-caching-location-update

## I had to do this step, i don't know if you will need to,try to run the estimation code without this first and if you get the error:
## "DEFUNCT: As of ExperimentHub (>1.17.2), default caching location has changed" then you need to run the following.
rappdirs::user_cache_dir(appname="ExperimentHub")
tools::R_user_dir("ExperimentHub", which="cache")

moveFiles<-function(package){
  olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
  newdir <- tools::R_user_dir(package, which="cache")
  dir.create(path=newdir, recursive=TRUE)
  files <- list.files(olddir, full.names =TRUE)
  moveres <- vapply(files,
                    FUN=function(fl){
                      filename = basename(fl)
                      newname = file.path(newdir, filename)
                      file.rename(fl, newname)
                    },
                    FUN.VALUE = logical(1))
  if(all(moveres)) unlink(olddir, recursive=TRUE)
}

package="ExperimentHub"
moveFiles(package)

## Now we can run the estimation code

Cell_type_EPIC<- estimateCellCounts2(rgset, compositeCellType = "Blood", 
                                     processMethod = "auto", cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), 
                                     referencePlatform = "IlluminaHumanMethylationEPIC")

## Now to map cell type as a variable to the targets file, first change to dataframe to prepare a merge
Cell_type_EPIC <- as.data.frame(Cell_type_EPIC)
Cell_type_EPIC$counts.total <-  rowSums (Cell_type_EPIC)
Cell_type_EPIC$GID <- rownames(Cell_type_EPIC)
rownames(Cell_type_EPIC) <- c()
View(Cell_type_EPIC)
save(Cell_type_EPIC,file="~/processing/R Data/Core/CellTypeEpicEstimation_04OCT2021.Rdata")

## Adding sample cell type to the targets file, THIS CAN STUFF UP THE ORDER OF THINGS SO DOUBLE CHECK EVERYTHING IS OK
targets <- merge(targets,Cell_type_EPIC, by="GID")
View(targets)
save(targets,file="~/processing/R Data/Core/targets_cell_04OCT2021.Rdata")
write.table(targets, file="~/processing/R Data/CSV Files/targetswithcell_04OCT2021.csv",sep=",",col.names=NA)





############################################ Principal Component Analysis or PCA ##################################################
## I have done some reading which says that you shouldn't do a PCA on categorical variables unless binary, but continuous ok.
## We have been ok publishing so far and this was from an older analysis pipeline that did the same thing but im sure one day a mathamatician will review and say this is wrong.
## Fact of the matter is i dont think the pipeline we had actually did a full PCA, it just looked at the PCA data, we adjust for obvious sources of variation in models anyway.

## An alternative would be to do a Multiple Factor Analysis.
## https://stats.stackexchange.com/questions/5774/can-principal-component-analysis-be-applied-to-datasets-containing-a-mix-of-cont

##Calculate principal component values, save as an object called loadings 

fit_PCA <- prcomp(t(m), center=TRUE, scale=TRUE, retx=TRUE)
loadings = fit_PCA$x


## To see a plot of the first 10 principal components and their contribution to the variation in a SCREE plot
options(scipen=999)
pdf(file="~/processing/Output/Scree PCA plot of 1st 10 principal components 04OCT2021.pdf")
plot(fit_PCA, type="lines", main="Scree PCA plot of 1st 10 principal components 04OCT2021")
dev.off()

## In my plot it looks like there are three sources of major variation.

## Now we build a principle component model


#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
#install.packages("WGCNA")
#library(WGCNA)
nGenes = nrow(m)
nSamples = ncol(m)
#well <- read.csv(file="~/wellID.csv",header=TRUE)
#targets <- merge(targets,well, by="GID")
#View(targets)

#chip <- read.csv(file="~/chip.csv",header=TRUE)
#targets <- merge(targets,chip, by="GID")
#View(targets)

#targets$Cur_PTSD_No <- ifelse(targets$CurrentPTSD=="No", "0",
#                      ifelse(targets$CurrentPTSD=="Yes", "1", NA))

#targets$Living_No <- ifelse(targets$Living=="City", "2",
#                              ifelse(targets$Living=="Village", "1", NA))

## Create a  data frame to add cell type and other  variables, this is just a simple model
datTraits = data.frame(Age=as.numeric(targets$Age_Child),
                       MoAge=as.numeric(targets$Age_Mother),
                       Cortisol=as.numeric(targets$Cortisol),
                       CD8T=as.numeric(targets$counts.CD8T),
                       CD4T=as.numeric(targets$counts.CD4T),
                       NK=as.numeric(targets$counts.NK),
                       Bcell=as.numeric(targets$counts.Bcell),
                       Mono=as.numeric(targets$counts.Mono),
                       Neu=as.numeric(targets$counts.Neu))



## Test correlation between principal components and data traits and find p values for correlations
options(scipen=3)
moduleTraitCor = cor(loadings[ ,c(1:20)], datTraits, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
loadings2 <- loadings[,c(1:20)]

png(file="~/processing/PCA/Data Heatmap of 20 principal components and various phenotypes 04OCT2021.png", height=720, width=720)
labeledHeatmap(Matrix = t(moduleTraitCor), xLabels = colnames(loadings2), yLabels= names(datTraits), xColorLabels=TRUE, yColorLabels=TRUE, colors=blueWhiteRed(6), textMatrix=t(textMatrix), setStdMargins=FALSE, cex.text=0.7, zlim=c(-1,1), main=paste("PCA-trait relationships")) 
dev.off()


## In my data, predicted cell type as well as maternal and offspring age have a great association with the differences in data. 
## These will need to be controlled for in further analysis (which we would anyway usually).





################################### Epigenetic Age estimation using Horvath's online clock calculator ################################### 

## You need a few things to run this generation, firstly the internet so you can make a log-in and can access this site https://dnamage.genetics.ucla.edu/new
## You can find information on the clock here: https://dnamage.genetics.ucla.edu/home
## Secondly you will need a couple of files i have made which will basically make a subset of 450 probes which are used in the generation of the clock
## I have conveniently uploaded them here for public use here: https://github.com/pfransquet/Epic-for-Horvath

setwd("~/processing/Clock")
## Download '450K_probes_missing_from_EPIC.Rdata', and 'HorvathClockProbes.Rdata' from https://github.com/pfransquet/Epic-for-Horvath and put them in your 'clock' folder
## Now we need to make a beta file with no probes removed, and normalise using the 'nOOb' method.
## Load in the rgset you saved, normalise (can take some time), map to genome and make a beta data set.
load("~/processing/R Data/Core/rgset04OCT2021.Rdata")
mset.noob <- preprocessNoob(rgset)
gmset.noob = mapToGenome(mset.noob)
beta.noob <- getBeta(gmset.noob)

Clockbetas <- as.data.frame(subset(beta.noob, rownames(beta.noob) %in% Clock450K_probe_list$ProbeID))
nrow(Clockbetas)

## Now need a separate targets file with ID, Age, Tissue and Sex included, which will be used as the "annotation file" for uploading to the clock.
## There is a list of tissue types that are predicted here https://dnamage.genetics.ucla.edu/sites/all/files/tutorials/TUTORIALonlineCalculator.pdf
## As I extracted DNA from whole blood, my Tissue type will be "Blood WB".
targets$Tissue<- "Blood WB"

## Find which rows you need and pick them out. In my case my ID's are row 2, Age is row 20, Sex is row 19 and Tissue is row 40
## Note the ID's column needs to be the same as the ID's you used when making the rgset so that it links up.
annotation_file <- targets[,c(2,20,19,40)]
write.csv(annotation_file, file = "annotation_file.csv",row.names = FALSE)


## Now we need to fix the file for upload. It is easier just to do this in excel than for me to write how to change it here as each data set will differ.
## The annotation file needs headings 'SampleID' 'Age' 'Tissue' and 'Female' with capitals at the start of each word, as if it is not formatted this way it will not work.
## "Female" (needs to have values 1 for female, 0 for male, NA for missing info)


## As there are a stack of probes from the 450k data set that aren't on the EPIC array, we have to put NAs for the missing probes which will be used for the clock.
## You can fill blanks as NA using the following code, combining the beta data with the missing probe data frame, using rbind.fill.
AllProbes <- rbind.fill(Clockbetas, missingprobe450K)
View(AllProbes)
## If you have a small data set (~80 samples or so) you can just continue on, save the "Allprobes" and this will be the "Methylation Data File" needed to upload to the clock website
write.csv(AllProbes, file = "betas.csv",row.names = FALSE)



### NOTE: I ORIGINALLY SPLIT MY DATA AS BELOW, BUT AS YOU ARE RUNNING DATA BASED ON THE WHOLE SET IT IS BETTER TO KEEP IT TOGETHER
### THE BETA's ACTUALLY 'ZIP" QUITE WELL TOGETHER, SO I JUST SUGGEST ZIPPING IT AND IF IT IS <800MB YOU ARE ALL GOOD!!


## However, the online clock calculator wants beta files to be <800MB so might need to split it up depending on size, in my case need to do in 4 batches (I have 260 samples) because once compile total file size is 1.4GB 
## Remember the probe ID is in the first place so need that in both files, it is tricky so just make sure the beta files line up with the annotation files you make below
## Again you might not need to do this, this is just and example of what i had to do.
betas1 <- AllProbes[,c(1:66)]
betas2 <- AllProbes[,c(1,67:132)]
betas3 <- AllProbes[,c(1,133:198)]
betas4 <- AllProbes[,c(1,199:262)]

## Also need to subset the targets file to match.
annotation_file1 <-annotation_file[c(1:65),]
annotation_file2 <-annotation_file[c(66:131),]
annotation_file3 <-annotation_file[c(132:197),]
annotation_file4 <-annotation_file[c(198:261),]

write.csv(betas1, file = "betas1.csv")
write.csv(betas2, file = "betas2.csv")
write.csv(betas3, file = "betas3.csv")
write.csv(betas4, file = "betas4.csv")
write.csv(annotation_file1, file = "annotation_file1.csv")
write.csv(annotation_file2, file = "annotation_file2.csv")
write.csv(annotation_file3, file = "annotation_file3.csv")
write.csv(annotation_file4, file = "annotation_file4.csv")

## I used excel now to change headings to work in the clock, including "ProbeID" for the beta sheets, 


## Now need to log in here https://dnamage.genetics.ucla.edu/new
## Clock wants beta file to be <800MB so might need to split it up depending on size, in my case need to do in a couple of batches
## If you have a big data set it can take a while to return, so to make sure your data is properly formatted
## You can just take the first 4 samples or so from the annotation (so sample rows including headings) and from beta (sample columns including probe ID), save them as "test" csv and upload those.
## IF you get a response then all good to go.

## Put in your email address
## Upload your betas in the "Methylation Data File" section
## Upload the annotation file you made in the "Sample Annotation File" section
## Select the tick boxes for 'Normalize data' and 'Advanced Analysis'
## Hit submit and stay on page until you get confirmation that everything has gone though, it can take a while depending on your internet upload speed.
## It can also crash and you have to do it again, i find if i am on the VPN i have troubles.

## Once successful the page will change and will show just this message:

## "File uploaded.
## You will receive results in the email when the processing is done.
## Go Back"

## If it has worked you will get results in an email. Even with big files if you do not get an email within a few hours, it has not worked and you will not get confirmation.
## This is why i test with smaller files to make sure all the formatting is correct.

## Once you get the email you will have the clock data, with over 100 variables including comments, cell and tissue estimation, as well as several clocks and other epigenetic measures.
## Look at the warnings and make sure there is nothing serious, sometimes you might have discordant sex, in my example some sample tissue type dont match up, but it is all labled as blood so ok.
## I would keep this raw file handy, then make another one with the variables that you need/want, there are some intersting things there, but in general if doing clock work on blood, i would keep:
## DNAmAge (which is the hovath measure), AgeAccelerationDiff,	AgeAccelerationResidual, SampleID.1 (you will need to link back to your targets file), 
## DNAmAgeHannum, DNAmPhenoAge	DNAmAgeSkinBloodClock, DNAmGrimAge, DNAmPhenoAgeAdjAge	DNAmAgeSkinBloodClockAdjAge, DNAmGrimAgeAdjAge,	IEAA,	EEAA,	IEAA.Hannum, AgeAccelerationResidualHannum,	AgeAccelGrim,	and AgeAccelPheno.
## More info on what each measure is can be found on the clock site.

## Change the 'SampleID.1' variable to whatever the ID variable name is in your targets file (for me it is GID), and then save this new file as "clockoutput.csv" and upload back into the "~/processing/Clock" folder.

## Now we want to link up the clock output back to the targets file using the ID variable.
Clock_Output <- read.csv(file="clockoutput.csv",header = TRUE)
targets <- merge(targets, Clock_Output, by="GID") ## put your identifier here where i have put 'GID'
## Just be aware that this will reorder your targets file by your identifier.
## Now just for records
write.csv(targets, file = "targets_with_clock.csv")
save(targets,file="~/processing/R Data/Core/targets_cell_clock_04OCT2021.Rdata")



############################################################ ANALYSIS USING LIMMA ############################################################ 
## Limma is used for linear models, using the 'm' data.
## First we do a case control, with no other variables in the model. I am using pregnancy PTSD in this example, a binary measure.
## I use '_raw' so i know this is the un adjusted linear model with no other variables.
Pregnancy_PTSD <- factor(targets$status)
design <- model.matrix(~Pregnancy_PTSD)
head(design)
fit <- lmFit(m,design)
fit2 <- eBayes(fit)
## Can get a summary to see if any probes are significant 
summary(decideTests(fit2))
## Now we make a table of the linear model results.
limmatable_raw = topTable(fit2, num=Inf, coef=2, adjust.method="BH")
View(limmatable_raw)
## Some information on the output can be found by typing
??limma::topTable
## and selecting the limma::topTable from the help box the pops up. You can also see here the other parameters you can enter when making the table.

## Now lets make linear models which adjust for other variables.
## I think here you should consider PCA or MFA output, and include in the model those variables that are shown to be highly associated with data variability above the outcome of interest or any known covariates.
## Here are two models one which adjusts for Age, Mothers age, chip well, batch, and sex, and another that includes those as well as all cell types.
## I put a VA on the end so i know it is a 'variable adjusted' model

## Just like to note which models are adjusting for what e.g. 
## VA1: Age, sex, Mothers age, chip well, batch.
## VA2: Age, sex, Mothers age, chip well, batch, and all cell types.

## Note: If you get an error "Error in lmFit(m, designcovar) : row dimension of design doesn't match column dimension of data object", it means that there might be data missing from one of the vars. 
## If any of your variables are missing data the linear model will not work, impute or remove the sample from both the m and targets datasets
## In my example one participant is missing mother data, it needs to be taken out of m, beta and targets "375_P"
## Just make sure to check that everything lines up as R is happy just to jumble everything randomly
which(colnames(m)=="375_P" )
m2 <- m[,-104]
which(colnames(beta)=="375_P" )
beta2 <- beta[,-104]
rownames(targets) <- targets$GID
which(rownames(targets)=="375_P")
targets2 <- targets[-104,]

## Model VA1 adjusting for Age, sex, Mothers age, chip well, and batch, i think binary files you dont need to put 'as.factor' but put it for other factors or the model will think it is continuous.
## Also note if you use this code, it will override some of the building blocks from the first model, so if you want bits and pieces you either have to re-run the first model code, or need to rename the data being made in this code
Pregnancy_PTSD <- factor(targets2$status)
designcovar <- model.matrix(~Pregnancy_PTSD + targets2$Age_Child + targets2$Gender_Child + targets2$Age_Mother + as.factor(targets2$Chip_well) + as.factor(targets2$Chip))
head(designcovar)
fitcovar <- lmFit(m2,designcovar)
fit2covar <- eBayes(fitcovar)
summary(decideTests(fit2covar))
## Pregnancy PTSD again was not sig for me, however chip 15 was, if you want data from a specific factor, you can change the 'coef=' number, e.g. the first factor of interest is the second factor in the model so thats why
## it is generally coef=2, if you want something else, for example there is one CpG associated with gender, and gender is my 4th factor, i can put 'coef=4' which give output specific to that factor.
limmatable_covarRun1 = topTable(fit2covar, num=Inf, coef=2)
## It should order the CpG's by adjusted p val, if not you can run the following code(of if you want to order by a different variable)
## limmatable_covarRun1 <- limmatable_covarRun1[order(limmatable_covarRun1$adj.P.Val), , drop = FALSE]
View(limmatable_covarRun1)

## Put everything from covar model 1 into a table including betas (so you can see % methylation), delta betas, epic annotation and limma stats
annEPICSub <- ann.EPIC[match(rownames(m2),ann.EPIC$Name), c(1:4,12:19,24:ncol(ann.EPIC))]
## Now subset based on data in the annEPICSub just made
DMP_VA1 = topTable(fit2covar, num=Inf, coef=2, genelist = annEPICSub)
View (DMP_VA1)

## Make deltabeta data (i.e the average difference between groups at each probe)
## I had to do this again because i had to subset the dataframes due to mothers missing age.
Control2=beta2[,targets2$status == "C"]
Case2=beta2[,targets2$status == "P"]
deltabeta2=rowMeans(Case2) - rowMeans(Control2)
dbeta2=data.frame(ID=names(deltabeta2), deltabeta2=deltabeta2)
View(dbeta2) ## I just like to compare the order of the deltabeta and the beta files to make sure they line up

allbeta=merge(beta2,dbeta2, by=0, all=TRUE)
rownames(allbeta) = allbeta$Row.names
allbeta$Row.names = NULL

## Now merge all data together, this will include the participant beta data, the annotation data, and the output from the model all together.
DMPbeta=merge(allbeta,DMP_VA1,by=0, all=TRUE)
## Now you can subset based on any condition, I will choose a stringent cutoff of unadjusted p value of 0.01
TopDMP_VA1 <- subset(DMPbeta, P.Value <= 0.01)
TopDMP_VA1 <-TopDMP_VA1[order(TopDMP_VA1$P.Value),] ## or you can use adj.P.Val to adjust by that variable.


## Now a second Model VA2 adjusting for Age, sex, Mothers age, chip well, batch, and all estimated cell types.
Pregnancy_PTSD <- factor(targets2$status)
designcovar <- model.matrix(~Pregnancy_PTSD + targets2$Age_Child + targets2$Gender_Child + targets2$Age_Mother + as.factor(targets2$Chip_well) + as.factor(targets2$Chip) 
                            + targets2$counts.Bcell + targets2$counts.CD4T  + targets2$counts.CD8T + targets2$counts.Mono + targets2$counts.Neu + targets2$counts.NK)
head(designcovar)
fitcovar <- lmFit(m2,designcovar)
fit2covar <- eBayes(fitcovar)
summary(decideTests(fit2covar))
limmatable_covarRun2 = topTable(fit2covar, num=Inf, coef=2)
View(limmatable_covarRun2)
annEPICSub <- ann.EPIC[match(rownames(m2),ann.EPIC$Name), c(1:4,12:19,24:ncol(ann.EPIC))]
DMP_VA2 = topTable(fit2covar, num=Inf, coef=2, genelist = annEPICSub)
View (DMP_VA2)
Control2=beta2[,targets2$status == "C"]
Case2=beta2[,targets2$status == "P"]
deltabeta2=rowMeans(Case2) - rowMeans(Control2)
dbeta2=data.frame(ID=names(deltabeta2), deltabeta2=deltabeta2)
View(dbeta2)

allbeta=merge(beta2,dbeta2, by=0, all=TRUE)
rownames(allbeta) = allbeta$Row.names
allbeta$Row.names = NULL
DMPbeta=merge(allbeta,DMP_VA2,by=0, all=TRUE)
TopDMP_VA2 <- subset(DMPbeta, P.Value <= 0.01)
TopDMP_VA2 <-TopDMP_VA2[order(TopDMP_VA2$P.Value),] ## or you can use adj.P.Val to adjust by that variable.
View(TopDMP_VA2)


############################# CATAGORICAL and CONTINUOUS MODELS #################################
## In your data set you may have a categorical, rather than binary, variable that you want to look at
## Here is an example using Low, normal and High cortisol levels (so a three catagory factor, in my targets file called "TertileCort")
## Note that if you have 3 categories, to correctly label based on what you want to compare to, in TertileCort, normal is coded a 0, low is 1 and 2 is high.
## The model will automatically compare to 0. So this model will make coef=2 normal vs low cortisol, and coef=3 normal vs high cortisol.
## The rest is much like the above analysis

##  Cortisol Tertile Model, with cell type
Child_Cort <- factor(targets2$TertileCort)
designcovar <- model.matrix(~Child_Cort + targets2$Age_Child + targets2$Gender_Child + targets2$Age_Mother + targets2$Chip_well + targets2$Chip + targets2$counts.CD8T  + targets2$counts.CD4T + targets2$counts.NK +
                              targets2$counts.Bcell + targets2$counts.Mono + targets2$counts.Neu)
head(designcovar)  
fitcovar <- lmFit(m2,designcovar)
fit2covar <- eBayes(fitcovar)
summary(decideTests(fit2covar))

## As it is catagorical, i need to split the output into normal vs low and normal vs high
## First normal vs low
limmatable_Low_Cort = topTable(fit2covar, num=Inf, coef=2)
limmatable_Low_Cort <- limmatable_Low_Cort[order(limmatable_Low_Cort$adj.P.Val), , drop = FALSE]
View(limmatable_Low_Cort)

annEPICSub <- ann.EPIC[match(rownames(m2),ann.EPIC$Name), c(1:4,12:19,24:ncol(ann.EPIC))]
DMP_lowcort = topTable(fit2covar, num=Inf, coef=2, genelist = annEPICSub)
View (DMP_lowcort)

## Now normal vs high
limmatable_High_Cort = topTable(fit2covar, num=Inf, coef=3)
limmatable_High_Cort <- limmatable_High_Cort[order(limmatable_High_Cort$adj.P.Val), , drop = FALSE]
View(limmatable_High_Cort)

DMP_highcort = topTable(fit2covar, num=Inf, coef=3, genelist = annEPICSub)
View (DMP_highcort)

## Delta betas are now going to be a little trickier as you cannot just summarise cases vs controls, we need averages of normal vs low, and normal vs high.
Normal=beta2[,targets2$TertileCort== "0"]
Low=beta2[,targets2$TertileCort== "1"]
High=beta2[,targets2$TertileCort== "2"]

dbnvl=rowMeans(Low) - rowMeans(Normal)
dbnvh=rowMeans(High) - rowMeans(Normal)
dbetalow=data.frame(ID=names(dbnvl), dbnvl=dbnvl)
dbetahigh=data.frame(ID=names(dbnvh), dbnvh=dbnvh)
View(dbetalow)
View(dbetahigh)

## Prepare the normal vs low dataset
lowallbeta=merge(beta2,dbetalow, by=0, all=TRUE)
rownames(lowallbeta) = lowallbeta$Row.names
lowallbeta$Row.names = NULL
DMPlowbeta=merge(lowallbeta,DMP_lowcort,by=0, all=TRUE)
TopDMP_lowcort <- subset(DMPlowbeta, P.Value <= 0.01)
TopDMP_lowcort <-TopDMP_lowcort[order(TopDMP_lowcort$P.Value),] ## or you can use adj.P.Val to adjust by that variable.
View(TopDMP_lowcort)
## Now to remove samples that are classed as high (as they were not used in comparisons)
highsamps <- subset(targets2,TertileCort=="2")
highsamps <- highsamps$GID
TopDMP_lowcort = TopDMP_lowcort[,!(names(TopDMP_lowcort) %in% highsamps)]


## Now the normal vs high dataset
highallbeta=merge(beta2,dbetahigh, by=0, all=TRUE)
rownames(highallbeta) = highallbeta$Row.names
highallbeta$Row.names = NULL
DMPhighbeta=merge(highallbeta,DMP_highcort,by=0, all=TRUE)
TopDMP_highcort <- subset(DMPhighbeta, P.Value <= 0.01)
TopDMP_highcort <-TopDMP_highcort[order(TopDMP_highcort$P.Value),] ## or you can use adj.P.Val to adjust by that variable.
View(TopDMP_highcort)
## Now to remove samples that are classed as low so you only have the normal and high
lowsamps <- subset(targets2,TertileCort=="1")
lowsamps <- lowsamps$GID
TopDMP_highcort = TopDMP_highcort[,!(names(TopDMP_highcort) %in% lowsamps)]



## Now lets do a continuous model, in my data i have cortisol measures which are continuous.
## I haven't had to really interpret this model properly so you may need to check your self.

## first lets do an unadjusted model
designcovar <- model.matrix(~targets2$Cortisol)
head(designcovar)  
fitcovar <- lmFit(m2,designcovar)
fit2covar <- eBayes(fitcovar)
summary(decideTests(fit2covar))
limmatable_Cort_cont = topTable(fit2covar, num=Inf, coef=2)
limmatable_Cort_cont <- limmatable_Cort_cont[order(limmatable_Cort_cont$adj.P.Val), , drop = FALSE]
View(limmatable_Cort_cont)

## If you want to add betas/cut the table down you can follow the instructions from above.
## Now lets do an adjusted model
designcovar <- model.matrix(~targets2$Cortisol +  targets2$Age_Child + targets2$Gender_Child + targets2$Age_Mother + targets2$Chip_well + targets2$Chip + targets2$counts.CD8T  + targets2$counts.CD4T + targets2$counts.NK +
                              targets2$counts.Bcell + targets2$counts.Mono + targets2$counts.Neu)

head(designcovar)  
fitcovar <- lmFit(m2,designcovar)
fit2covar <- eBayes(fitcovar)
summary(decideTests(fit2covar))
limmatable_Cort_cont_adj = topTable(fit2covar, num=Inf, coef=2)
limmatable_Cort_cont_adj <- limmatable_Cort_cont_adj[order(limmatable_Cort_cont_adj$adj.P.Val), , drop = FALSE]
View(limmatable_Cort_cont_adj)

## You can subset and make a couple of different data frames based on cutoff that are more likely to be meaningful, starting with Pvalue, then choosing a log fold change cutoff
limmatable_Cort_cont_adj <- subset(limmatable_Cort_cont_adj, limmatable_Cort_cont_adj$P.Value <= 0.01)

limmatable_Cort_cont_adj_up <- subset(limmatable_Cort_cont_adj, logFC > 0.0005)
limmatable_Cort_cont_adj_down <- subset(limmatable_Cort_cont_adj, logFC < -0.0005)

limmatable_Cort_cont_adj_up <- limmatable_Cort_cont_adj_up[order(limmatable_Cort_cont_adj_up$logFC,decreasing = TRUE), , drop = FALSE]
limmatable_Cort_cont_adj_down <- limmatable_Cort_cont_adj_down[order(limmatable_Cort_cont_adj_down$logFC), , drop = FALSE]

Top_DMP_Cort_cont_adj_by_FC <-rbind(limmatable_Cort_cont_adj_up, limmatable_Cort_cont_adj_down)


############################ After analysis data visualisation ##########################
## Sometimes it is handy to look at the data, particularly the high sig ones (or large effect size) to see what the data looks like
## Often you might get patterns that suggest genetic variation that is driving the methylation signal.
## I like using scatter plots, lets use the output from the model we just ran, keep in mind this will plot raw data, and not adjusted data.
## This will make a list of the pvaues (as ordered above) and then make plots based on the top 200 pvals, and i will use the negative fold change.
limmatable_Cort_cont_adj_down$cpg = row.names(limmatable_Cort_cont_adj_down)
limmatable_Cort_cont_adj_down <- limmatable_Cort_cont_adj_down[order(limmatable_Cort_cont_adj_down$P.Value), , drop = FALSE]
list <- limmatable_Cort_cont_adj_down[c(1:200), c(7)]
pdf(file="~/processing/Top 200 DMPs ranked by -FC Scatter Box.pdf", width = 14)
par(mfrow=c(2,5))
plotCpg(dat=beta2,pheno=targets2$Cortisol,type="continuous", fitLine = TRUE, cpg=list, ylab = "Beta Values")
dev.off()
## You can make a csv of these cpgs if you want too 
## write.csv(list, file = "~/processing/List of top 200 DMPs ranked by -FC.csv")


## And again for the top 200 pvals but using positive correlation 
limmatable_Cort_cont_adj_up$cpg = row.names(limmatable_Cort_cont_adj_up)
limmatable_Cort_cont_adj_up <- limmatable_Cort_cont_adj_up[order(limmatable_Cort_cont_adj_up$P.Value), , drop = FALSE]
list <- limmatable_Cort_cont_adj_up[c(1:200), c(7)]
pdf(file="~/processing/Top 200 DMPs ranked by posFC Scatter Box.pdf", width = 14)
par(mfrow=c(2,5))
plotCpg(dat=beta2,pheno=targets2$Cortisol,type="continuous", fitLine = TRUE, cpg=list, ylab = "Beta Values")
dev.off()

## If you want to just do one probe, for example the highest ranked cpg with adjusted p val was cg24650748 in the negative correlations
pdf(file="~/processing/Top probe by adjusted pval -FC.pdf", width = 7)
par(mfrow=c(1,1))
plotCpg(dat=beta2,pheno=targets2$Cortisol,type="continuous", fitLine = TRUE, cpg="cg24650748", ylab = "Beta Values")
dev.off()


## Now lets have a look at some of the significant probes (top 200) from the case control type analysis (Pregnancy PTSD). Just make sure you are using the right limma table, you might have to run the model again.
limmatable_covarRun2$cpg = row.names(limmatable_covarRun2)
limmatable_covarRun2 <- limmatable_covarRun2[order(limmatable_covarRun2$P.Value), , drop = FALSE]
list <- limmatable_covarRun2[c(1:200), c(7)]
pdf(file="~/processing/Top 200 DMPs, pregnancy PTSD (adj for Age, sex, Mothers age, chip well, batch, and all estimated cell types).pdf", width = 14)
par(mfrow=c(2,5))
plotCpg(dat=beta2,pheno=targets2$Pregnancy_PTSD,type="categorical", cpg=list, ylab = "Beta Values")
dev.off()

## If you want to just do one probe...
pdf(file="~/processing/Top probe, pregnancy PTSD (adj for Age, sex, Mothers age, chip well, batch, and all estimated cell types).pdf", width = 7)
par(mfrow=c(1,1))
plotCpg(dat=beta2,pheno=targets2$Pregnancy_PTSD,type="categorical", cpg="cg14404316", ylab = "Beta Values")
dev.off()


## There are fancier ways to plot CpGs too using ggplot, they take a bit of data manipulation though.
## There are some examples below as part of the cate method.

#######################################################  MDS plots based on methylation #####################################################
## We will use the limma table to make some examples here, make sure your limma table is in order of best probes based on P Value
limmatable_covarRun1 <- limmatable_covarRun1[order(limmatable_covarRun1$P.Value), , drop = FALSE]
limmatable_covarRun1$probe <- rownames(limmatable_covarRun1)

## MDS plot for top 20 probes
pdf(file="~/processing/DMPs/Top 20 DMPs MDS.pdf", height = 10, width = 10)
par(mfrow=c(1,1))
top20list <- limmatable_covarRun1[c(1:20), c(7)]
top20beta<-beta[top20list,]
mdsPlot(top20beta,numPositions=nrow(top20beta),sampNames=NULL,pch=19,sampGroup=targets$Pregnancy_PTSD)
dev.off()

## Can do the same for any amount of probes
## MDS plot for top 50 probes
pdf(file="~/processing/DMPs/Top 50 DMPs MDS.pdf", height = 10, width = 10)
par(mfrow=c(1,1))
top50list <- limmatable_covarRun1[c(1:50), c(7)]
top50beta<-beta[top50list,]
mdsPlot(top50beta,numPositions=nrow(top50beta),sampNames=NULL,pch=19,sampGroup=targets$Pregnancy_PTSD)
dev.off()

## Or even all sig results
## All sig 
pdf(file="~/processing/DMPs/All Sig DMPs MDS.pdf", height = 10, width = 10)
par(mfrow=c(1,1))
top_sig <- subset(limmatable_covarRun1, P.Value <= 0.05)
top_sig_list <- top_sig[, c(7)]
top_sig_beta<-beta[top_sig_list,]
mdsPlot(top_sig_beta,numPositions=nrow(top_sig_beta),sampNames=NULL,pch=19,sampGroup=targets$Pregnancy_PTSD)
dev.off()



## You could also manupulate the above code to take your Top DMP data, sort by negative or positive deltabeta to do similar plots with a focus on effect size, or just positive/negative effect sizes etc.



############################### Alternative analysis method 'cate' ##################################
## https://cran.r-project.org/web/packages/cate/index.html
## https://cran.r-project.org/web/packages/cate/cate.pdf
## Because EWAS datasets are so vast there is lots of possibility for noise, and with small data sets you might not get any probes which pass adjustement for multiple testing.
## This method 'cate' removes unwanted variation which increases 'real' signals and perhaps gives better chance of getting something significant amoungst the noise.
## This is what i used to do the Cortisol analysis instead of the limma models above.
## You need to subset M  based on whatever model you are using, so if you are using all M data then no need to subset.

## I made categorical variables binary where i can so it isn't pulling the data around too much in adjustment
## Here is an example to make categories Marital and Education into binary
targets2$Married <- ifelse(targets2$Marital=="1", "Married", "Other")
targets2$SecEdu <-ifelse(targets2$Education=="3"|targets2$Education=="4", "secondary_up", "primary_down")

## This will subset based on the sample IDs
samples <- targets2$GID
m_cate <- m[,samples]
beta_cate <- beta[,samples]
## cate' wants the outcome (i.e. methylation) data in participant rows and methylation columns rather than the usual participant columns and methylation rows, so it first needs to be transposed. 
## Again we are using m values over beta's as they are reported to be better than beta when carrying out stats.
## Also the number of participants has to be the same for both data sets
mt <-t(m_cate)

## Check that there are the same number of rows in both
nrow(mt)
nrow(targets2)
## If they aren't, make sure your row names are the same identifier as the the rownames in the m file made 
rownames(mt)
rownames(targets2)
## I summarised the column names here so i know what might be of interest without having to bring up the targets table every time
colnames(targets2)

## Cell estimates
#"counts.CD8T"
#"counts.CD4T"
#"counts.NK"
#"counts.Bcell"
#"counts.Mono"    
#"counts.Neu"

## As well as technical variables
#"Chip_well"      
#"Chip"           

##And participant characteristic variables, 
#"Gender_Child"   
#"Age_Child"      
#"Age_Mother"    
#"Living"
#"Education"
#"Marital"  
#"Smoke_Preg"     
#"ID" 

##As well as the binary ones we made
#"Married"
#"SecEdu"

## And the variables of interest to this example
#"Pregnancy_PTSD" (a binary case/control)
#"Cortisol" (Continuous variable)
#"TertileCort" (Catagorical variable)

## Basically with cate you have to make the targets file your model, so whatever variables there are in there it will use.
## You can make a separate targets file having only the variables you might be interested in, so you do not have to keep removing variables you will never need.
## So in this case i will remove all data in the targets file i have no need for and call it targets_cate, which will be my foundation target file for these models.
## Go in and write down the column number of the values you do not want and add them to the code to remove variables while making a new data set.
targets_cate <- targets2[-c(1,2,9:12,15:18,22,23,28,29,33:43)]
View(targets_cate)


## For the first example, lets use a similar model as the limma model 'VA2', a prenatal PTSD case control, adjusting for Age, sex, Mothers age, and all cell types.
## Like the limma models, if any participant is missing data for a variable the model will just error.
## You should specify which variables are categorical factors or it will think it is continuous, if they are binary it will be ok as is. 
## If you are adjusting for chip/well, you have to specify in the model it is catagorical, so i think you have to also list all the other variables, which we will do below in a second example.
## First have to remove any variables that are not need in the model (keeping only prenatal PTSD, Age, sex, Mothers age, and all cell types)
targets_VA2 <- targets_cate[-c(7,8,13:20)]
View(targets_VA2)

## This next line can take some time (1 hr per 100 samples or so), you can make it run faster by chaning nRepeat lower, but the model will be 'less accurate'
## This is basically telling the program to use pregnancy ptsd as the factor of interst, and all other variables in the targets table as possible confounding
factor.num <- est.confounder.num(~ Pregnancy_PTSD |. - Pregnancy_PTSD + 0 , targets_VA2, mt,  method = "bcv", bcv.plot = TRUE, rmax = 30, nRepeat = 100)
## See the number of estimated latent factors to adjust
factor.num$r
## now get the results of the model, this only takes a couple of minutes
cate.results <- cate(~ Pregnancy_PTSD | . - Pregnancy_PTSD +0  , targets_VA2, mt, r = factor.num$r)
cate.results$alpha.p.value
## now check which results passed a cut off of your liking
which(p.adjust(cate.results$beta.p.value, "BH") < 0.1) 
## For this model, nothing passed at less conservative <0.1
## In the following models there are some which i will use for an example on what to do next.

## Now lets run the same model, but with chip well and batch, which we have to list as categorical variables.
targets_VA3 <- targets_cate[-c(13:20)]
View(targets_VA3)
factor.num <- est.confounder.num(~ Pregnancy_PTSD | counts.CD8T + counts.CD4T + counts.NK + counts.Bcell + counts.Mono + counts.Neu + as.factor(Chip_well) + as.factor(Chip) 
                                 + as.factor(Gender_Child) + Age_Child + Age_Mother + 0 , targets_VA3, mt,  method = "bcv", bcv.plot = TRUE, rmax = 30, nRepeat = 100)
factor.num$r
## In this model there is apparently no latent factors to adjust for, i include in model anyway but probably no need
cate.results <- cate(~ Pregnancy_PTSD | counts.CD8T + counts.CD4T + counts.NK + counts.Bcell + counts.Mono + counts.Neu + as.factor(Chip_well) + as.factor(Chip) 
                     + as.factor(Gender_Child) + Age_Child + Age_Mother + 0 , targets_VA3, mt, r = factor.num$r)
cate.results$alpha.p.value
## now check which results passed a cut off of your liking
which(p.adjust(cate.results$beta.p.value, "BH") < 0.1) 
## This time i got one probe, ill show you how to use this data further below


####################################################### CORT CONTINUOUS and CORT Catagorical MODELS ###################################################################
##Now lets do a continuous Model with cell, ages, sex, and demographic data, but i leave PREG PTSD out because i wanted to stratify after 
targets_cort <- targets_cate[-c(7,11,14,15,17)]
factor.num <- est.confounder.num(~ Cortisol |. - Cortisol + 0 , targets_cort, mt,  method = "bcv", bcv.plot = TRUE, rmax = 30, nRepeat = 100)
## Just for an example, when i ran this i got the error:
#  Error in cbind(X.nuis, X.primary) : 
#  number of rows of matrices must match (see arg 2)
## This is basically suggesting there is data missing from somewhere, in my case two participants are missing marital status and one smoking status. 
## If you want to adjust for these you must impute the data or remove the sample from the data. I am just going to remove the variables so the model runs
targets_cort <- targets_cate[-c(7,8,11,13:15,17:20)]
## Try again
factor.num <- est.confounder.num(~ Cortisol |. - Cortisol + 0 , targets_cort, mt,  method = "bcv", bcv.plot = TRUE, rmax = 30, nRepeat = 100)
factor.num$r
cate.results <- cate(~ Cortisol | . - Cortisol +0  , targets_cort, mt, r = factor.num$r)
cate.results$alpha.p.value
which(p.adjust(cate.results$beta.p.value, "BH") < 0.1) 
## Here i just get one that passes the cut off.

##Same model but including sample well and batch as categorical variables
targets_cort <- targets_cate[-c(11,13:15,17:20)]
factor.num <- est.confounder.num(~ Cortisol | counts.CD8T + counts.CD4T + counts.NK + counts.Bcell + counts.Mono + counts.Neu + as.factor(Chip_well) + as.factor(Chip) 
                                 + as.factor(Gender_Child) + Age_Child + Age_Mother + 0 , targets_cort, mt,  method = "bcv", bcv.plot = TRUE, rmax = 30, nRepeat = 100)
factor.num$r
cate.results <- cate(~ Cortisol | counts.CD8T + counts.CD4T + counts.NK + counts.Bcell + counts.Mono + counts.Neu + as.factor(Chip_well) + as.factor(Chip) 
                     + as.factor(Gender_Child) + Age_Child + Age_Mother + 0 , targets_cort, mt, r = factor.num$r)
cate.results$alpha.p.value
which(p.adjust(cate.results$beta.p.value, "BH") < 0.15) 
## I got five <0.15 so i will go with that cut off for the example below
##[1]  46107 246187 404797 486172 508555
## Can also check different adjustment methods depending on what yout want to report, e.g.
which(p.adjust(cate.results$beta.p.value, "bonferroni") < 0.15)
## of which bonferroni give one [1] 486172

## Ok lets have a look at these probes! I've named the data frame res for 'Results'
## Basically this code manipulates the cate.results output for our use
res <- as.data.frame (cate.results$beta.p.value)
names(res)[names(res) == "Cortisol"] <- "pval"
res$beta <- cate.results$beta
res$t <- cate.results$beta.t
res$probe <- row.names(res)
row.names(res) <- NULL

sigprobe_rowlist <- as.character(which(p.adjust(cate.results$beta.p.value, "BH") < 0.15))
sigprobes <- as.character(res[sigprobe_rowlist,"probe"])
row.names(res) <- res$probe
View(res)
sigprobesdata <- res[sigprobes, ]
View(sigprobesdata) ## Now we can see which probes were significant, thier pval (non adjusted) as well as some other values.

## Now make a dataframe with the adjusted pvalues
adjustment <- as.data.frame(p.adjust(cate.results$beta.p.value, "BH"))
adjustmentNone <- as.data.frame(p.adjust(cate.results$beta.p.value, "none"))
adjustmentbonferroni <- as.data.frame(p.adjust(cate.results$beta.p.value, "bonferroni"))
adjustments<-cbind(adjustment,adjustmentbonferroni,adjustmentNone)
names(adjustments)[names(adjustments) == "p.adjust(cate.results$beta.p.value, \"BH\")"] <- "BH"
names(adjustments)[names(adjustments) == "p.adjust(cate.results$beta.p.value, \"none\")"] <- "Unadjusted"
names(adjustments)[names(adjustments) == "p.adjust(cate.results$beta.p.value, \"bonferroni\")"] <- "Bonferroni"
adjustments$probe <- row.names(res)

adjustments2 <- adjustments ## This will be all the probes in the model with the adjusted values
## As adjustments has numbered rows, can use sigprobe_rowlist, if has row lables which are cpgs then can use sig probes  ##adjustments2 <- adjustments2[sigprobes, ] ##
adjustments2 <- adjustments2[sigprobe_rowlist, ]
## Now you have a list of the sig probes with your chosen cut off.

rownames (adjustments2) <- c(adjustments2$probe)
Probedata <- merge(adjustments2, sigprobesdata, by=0, all=TRUE)
View(Probedata) ## here 'unadjusted' and 'pval' should be the same, and 'probe x', 'probe y', and 'row.names' should be the same.
## Now clean up after checking
names(Probedata)[names(Probedata) == "Row.names"] <- "Probe"
Probedata$probe.x = NULL
Probedata$probe.y = NULL
Probedata$pval = NULL
rownames(Probedata) <- c(Probedata$Probe) 
## Now merge with probe information from annotated EPIC file 'ann.EPIC', similar to the limma tables above.
probeinfo <- as.data.frame (ann.EPIC[sigprobes, ])
#row.names(probeinfo)  <- c(probeinfo$Name) # Incase they just come up as numbered rows
probebeta <- as.data.frame (beta_cate[sigprobes, ])
Probedata2 <- merge(Probedata, probeinfo, by="row.names", all=TRUE)
row.names(Probedata2)  <- c(Probedata2$Row.names)
Probedata2$Row.names = NULL
Probedata2 <- merge(Probedata2, probebeta, by="row.names", all=TRUE)
View(Probedata2)

## Now name it something more useful and save it
Continuous_cortisol_probes_adjusted_model <- Probedata2
write.csv(Continuous_cortisol_probes_adjusted_model, file = "Continuous_cortisol_probes_adjusted_model.csv")


## Now lets do a categorical model using tertile cortisol measures (Low, normal and High), adjusting for cell-type, batch, ages, and sex.
targets_cortcat <- targets_cate[-c(7,13:16,18:20)]
## Here is some code if you want to rename catagories, which might be useful for making graphs later
#targets_cortcat$cat <- ifelse(targets_cortcat$TertileCort=="0", "Normal", ifelse(targets_cortcat$TertileCort=="1", "Low", ifelse(targets_cortcat$TertileCort=="2", "High", NA))) 
factor.num <- est.confounder.num(~ as.factor(TertileCort) | counts.CD8T + counts.CD4T + counts.NK + counts.Bcell + counts.Mono + counts.Neu + as.factor(Chip) + as.factor(Gender_Child) 
                                 + Age_Child + Age_Mother + 0 , targets_cortcat, mt,  method = "bcv", bcv.plot = TRUE, rmax = 30, nRepeat = 100)
factor.num$r
cate.results <- cate(~ as.factor(TertileCort) | counts.CD8T + counts.CD4T + counts.NK + counts.Bcell + counts.Mono + counts.Neu + as.factor(Chip) + as.factor(Gender_Child) 
                     + Age_Child + Age_Mother + 0 , targets_cortcat, mt, r = factor.num$r)
cate.results$alpha.p.value
which(p.adjust(cate.results$beta.p.value, "BH") < 0.20) 
## [1]   16325  214297  304222  579643  749917  902947  955934 1196762 1196763 1279397 1291746 1503678
## This time i have 12 probes as i raised the bar, also because it has technically done two analyses there is an issue with numbering which I will explain below

## Just like the limma model, because there are several categories it splits the analysis so we have to capture both
res <- as.data.frame (cate.results$beta.p.value)
names(res)[names(res) == "as.factor(TertileCort)1"] <- "Norm_v_Low_pval"
names(res)[names(res) == "as.factor(TertileCort)2"] <- "Norm_v_High_pval"
res$beta <- cate.results$beta
res$t <- cate.results$beta.t
res$probe <- row.names(res)
row.names(res) <- NULL

sigprobe_rowlist <- as.character(which(p.adjust(cate.results$beta.p.value, "BH") < 0.20))
## Because there are two rows (norm vs low and norm vs high) the program just keeps counting from one row to another, resulting in the row identifiers being higher than the amount of rows.
## Note the number of rows in your data frame
nrow(res)
#760590
## to get the actual row you need (where the probe of interest is), you need to to subtract the number of rows in the data frame (so in my case 760590) from the number given if it is higher than the 
## number of rows in the data frame, so for example, probes at row numbers "16325"  "214297"  "304222"  "579643" and "749917" are ok, but for the others i need to subtract
## e.g, 902947-760590 = 142357
## 955934-760590 = 195344
## 1196762-760590 = 436172
## 1196763-760590 = 436173
## 1279397-760590 = 518807
## 1291746-760590 = 531156
## 1503678-760590 = 743088

## I went back to manually check, and it looks like it works out fine, also the ones >760590 (for me) will be from the second column, i.e. norm vs high
## you can make a seperate list with the propper numbers e.g.
sigprobe_no <-as.numeric(sigprobe_rowlist)
sigprobe_no <- sigprobe_no - nrow(res)
sigprobe_no ## ignore the negative '-' numbers and replace with the proper number from first row (numbers below 760590)

sigprobe_rowlist2 <- as.character(c("16325", "214297", "304222", "579643",  "749917",  "142357",  "195344",  "436172",  "436173",  "518807", "531156",  "743088"))
sigprobes <- as.character(res[sigprobe_rowlist2,"probe"])
row.names(res) <- res$probe
View(res)
sigprobesdata <- res[sigprobes, ]
View(sigprobesdata)

adjustment <- as.data.frame(p.adjust(cate.results$beta.p.value, "BH"))
adjustmentNone <- as.data.frame(p.adjust(cate.results$beta.p.value, "none"))
adjustmentbonferroni <- as.data.frame(p.adjust(cate.results$beta.p.value, "bonferroni"))
adjustments<-cbind(adjustment,adjustmentbonferroni,adjustmentNone)
names(adjustments)[names(adjustments) == "p.adjust(cate.results$beta.p.value, \"BH\")"] <- "BH"
names(adjustments)[names(adjustments) == "p.adjust(cate.results$beta.p.value, \"none\")"] <- "Unadjusted"
names(adjustments)[names(adjustments) == "p.adjust(cate.results$beta.p.value, \"bonferroni\")"] <- "Bonferroni"
adjustments$probe <- row.names(res)

adjustments2 <- adjustments
## As adjustments has numbered rows, can use sigprobe_rowlist, if has row lables which are cpgs then can use sig probes  ##adjustments2 <- adjustments2[sigprobes, ] ##
adjustments2 <- adjustments2[sigprobe_rowlist, ]
rownames (adjustments2) <- c(adjustments2$probe)

## Merge the probe data
Probedata <- merge(adjustments2, sigprobesdata, by=0, all=TRUE)
View(Probedata)
## Again here 'unadjusted' and 'pval' should be the same, and 'probe x', 'probe y', and 'row.names' should be the same.
## Now clean up after checking
names(Probedata)[names(Probedata) == "Row.names"] <- "Probe"
Probedata$probe.x = NULL
Probedata$probe.y = NULL
Probedata$pval = NULL
rownames(Probedata) <- c(Probedata$Probe)

## Now merge with probe information from annotated EPIC file 'ann.EPIC'

#probe_list <- c(Probedata$Probe)
probeinfo <- as.data.frame (ann.EPIC[sigprobes, ])
#row.names(probeinfo)  <- c(probeinfo$Name) # Incase they just come up as numbered rows
probebeta <- as.data.frame (beta_cate[sigprobes, ])
Probedata2 <- merge(Probedata, probeinfo, by="row.names", all=TRUE)
row.names(Probedata2)  <- c(Probedata2$Row.names)
Probedata2$Row.names = NULL
Probedata2 <- merge(Probedata2, probebeta, by="row.names", all=TRUE)
View(Probedata2)

## Now name it something more useful and save it
Catagorical_cortisol_probes_adjusted_model <- Probedata2
write.csv(Catagorical_cortisol_probes_adjusted_model, file = "Catagorical_cortisol_probes_adjusted_model.csv")

## Now lets do some visualization

install.packages("forcats")
library(forcats)
install.packages("ggpubr")
library(ggpubr)
## I'll take the top 4 sig probes (<0.15) from my continuous cortisol model for this example 
sigprobes
## So i will use
## "cg01975957" "cg07100000" "cg11932447" "cg17342241"
## Remember if you didn't name it something different from the code, the models will overwrite themselves so make sure you are using the right data!
## I'm sure there is a way to automate this but here is how to do it one by one.

## Get the beta values for the probes of interest
## This wont work further on if you haven't subset the beta file to the same participants as used in the cate models
## as i removed one sample, it also has to be removed from the beta file
## I have done this above at one stage but if you havent done it then you can do it like this
#samples <- rownames(targets_cate)
#beta_cate <- beta[,samples]
cg01975957 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg01975957"))
cg01975957 <- t(cg01975957)
cg07100000 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg07100000"))
cg07100000 <- t(cg07100000)
cg11932447 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg11932447"))
cg11932447 <- t(cg11932447)
cg17342241 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg17342241"))
cg17342241 <- t(cg17342241)

## Add them as a variable to the targets cate file (or make a new one if you want, as you might still need the targets_cate file as is)
targets_cate$cg01975957 <- cg01975957[,c(1)]
targets_cate$cg07100000 <- cg07100000[,c(1)]
targets_cate$cg11932447 <- cg11932447[,c(1)]
targets_cate$cg17342241 <- cg17342241[,c(1)]

## Make it out of 100 (so % rather than between 0 and 1)
targets_cate$cg01975957 <- targets_cate$cg01975957*100
targets_cate$cg07100000 <- targets_cate$cg07100000*100
targets_cate$cg11932447 <- targets_cate$cg11932447*100
targets_cate$cg17342241 <- targets_cate$cg17342241*100


## You can change any of the setting here
## I got to this after lots of tweeking but you might for example want bigger dots, difference coloured line, etc.

cg01975957_plot <-ggplot(targets_cate, aes(x=Cortisol, y=cg01975957)) + geom_point(size=1, color="black", shape=16, ) + geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + labs(x = "Cortisol (nmol/L)")
cg07100000_plot <-ggplot(targets_cate, aes(x=Cortisol, y=cg07100000)) + geom_point(size=1, color="black", shape=16, ) + geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + labs(x = "Cortisol (nmol/L)")
cg11932447_plot <-ggplot(targets_cate, aes(x=Cortisol, y=cg11932447)) + geom_point(size=1, color="black", shape=16, ) + geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + labs(x = "Cortisol (nmol/L)")
cg17342241_plot <-ggplot(targets_cate, aes(x=Cortisol, y=cg17342241)) + geom_point(size=1, color="black", shape=16, ) + geom_smooth(method=lm, se=FALSE, linetype="dashed", color="black") + labs(x = "Cortisol (nmol/L)")

ggarrange(cg01975957_plot, cg07100000_plot, cg11932447_plot, cg17342241_plot, labels = c("a)", "b)", "c)", "d)"), ncol = 4, nrow = 1)

## As it is just sitting in the plots, you can hit export and then choose how you want to export it.
## An easy way to keep quality is PDF, and you can move the user interface around to get it looking ok and then export it "Device Size"
## Once you have saved then you can remove it
dev.off()

## Now some visulisation for catagories of cortisol
## I had 12 probes associated with high or low cortisol, but i will just pick 4 at random to show
sigprobes
## I'll go with
#cg10095352
#cg25649709
#cg23987897
#cg23917477
targets_cate$cat <- ifelse(targets_cate$TertileCort=="0", "Normal", ifelse(targets_cate$TertileCort=="1", "Low", ifelse(targets_cate$TertileCort=="2", "High", NA)))  

cg10095352 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg10095352"))
cg10095352 <- t(cg10095352)
cg25649709 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg25649709"))
cg25649709 <- t(cg25649709)
cg23987897 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg23987897"))
cg23987897 <- t(cg23987897)
cg23917477 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg23917477"))
cg23917477 <- t(cg23917477)

targets_cate$cg10095352 <- cg10095352[,c(1)]
targets_cate$cg25649709 <- cg25649709[,c(1)]
targets_cate$cg23987897 <- cg23987897[,c(1)]
targets_cate$cg23917477 <- cg23917477[,c(1)]

targets_cate$cg10095352 <- targets_cate$cg10095352*100
targets_cate$cg25649709 <- targets_cate$cg25649709*100
targets_cate$cg23987897 <- targets_cate$cg23987897*100
targets_cate$cg23917477 <- targets_cate$cg23917477*100

cg10095352_plot <- ggplot(targets_cate, aes(x=as.factor(cat), y=cg10095352)) + geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1)+ geom_jitter(size=1, shape=16, width=0.075)+ scale_x_discrete(limits = c("Low", "Normal", "High"), labs(x = "Tertile Cortisol"))
cg25649709_plot <- ggplot(targets_cate, aes(x=as.factor(cat), y=cg25649709)) + geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1)+ geom_jitter(size=1, shape=16, width=0.075)+ scale_x_discrete(limits = c("Low", "Normal", "High"), labs(x = "Tertile Cortisol"))
cg23987897_plot <- ggplot(targets_cate, aes(x=as.factor(cat), y=cg23987897)) + geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1)+ geom_jitter(size=1, shape=16, width=0.075)+ scale_x_discrete(limits = c("Low", "Normal", "High"), labs(x = "Tertile Cortisol"))
cg23917477_plot <- ggplot(targets_cate, aes(x=as.factor(cat), y=cg23917477)) + geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1)+ geom_jitter(size=1, shape=16, width=0.075)+ scale_x_discrete(limits = c("Low", "Normal", "High"), labs(x = "Tertile Cortisol"))

## You can make a seperate plot for these catagorial model data
ggarrange(cg10095352_plot, cg25649709_plot, cg23987897_plot, cg23917477_plot, labels = c("e)", "f)", "g)", "h)"), ncol = 4, nrow = 1)


## As we made data plots for the above linar models, we can stuck it all together (very handy for making graphs for publications)
## More info found here https://rpkgs.datanovia.com/ggpubr/reference/ggarrange.html
ggarrange(cg01975957_plot, cg07100000_plot, cg11932447_plot, cg17342241_plot, cg10095352_plot, cg25649709_plot, cg23987897_plot, cg23917477_plot, labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"), ncol = 4, nrow = 2, hjust = -0.05)



## Here is an example of case control, looking at one CPG with seperate trendlines for case vs control
cg27193031 <- as.data.frame(subset(beta_cate, rownames(beta_cate)== "cg27193031"))
cg27193031 <- t(cg27193031)
targets_cate$cg27193031 <- cg27193031[,c(1)]
targets_cate$cg27193031 <- targets_cate$cg27193031*100
targets_cate$status <- ifelse(targets_cate$Pregnancy_PTSD=="0", "No", ifelse(targets_cate$Pregnancy_PTSD=="1", "Yes", NA))
cg27193031_plot1 <-ggplot(targets_cate, aes(x=Cortisol, y=cg27193031,  colour = status, linetype = status)) + geom_point(size=1, color="black", shape=16) + geom_smooth(method=lm, se=FALSE, color="black") + labs(x = "Cortisol (nmol/L)", y="BDNF (cg27193031)") +  scale_linetype_discrete(name="Prenatal\nPTSD", breaks=c("No", "Yes"), labels=c("No", "Yes")) + theme(legend.title = element_text(size = 6), legend.text=element_text(size=6), legend.key.size = unit(.5, "cm"))
cg27193031_plot1

##Again you can play around with the settings e.g. make the plot, then shrink the legend
cg27193031_plot2 <-ggplot(targets_cate, aes(x=Cortisol, y=cg27193031,  colour = status, linetype = status)) + geom_point(size=1, color="black", shape=16) + geom_smooth(method=lm, se=FALSE, color="black") + labs(x = "Cortisol (nmol/L)", y="BDNF (cg27193031)") +  scale_linetype_discrete(name="Maternal\nPTSD", breaks=c("No", "Yes"), labels=c("No", "Yes"))
cg27193031_plot2 + theme(legend.title = element_text(size = 6), legend.text=element_text(size=6), legend.key.width = unit(0.1,"cm"), legend.key.size = unit(.6, "cm"))

## Or take the legend away
cg27193031_plot <-ggplot(targets_cate, aes(x=Cortisol, y=cg27193031, colour = status, linetype = status)) + 
  geom_point(size=1, color="black", shape=16) + 
  geom_smooth(method=lm, se=FALSE, color="black") + 
  labs(x = "Cortisol (nmol/L)", y="BDNF (cg27193031)") + 
  scale_linetype_discrete(name="status", breaks=c("Control", "PTSD"), labels=c("Control", "PTSD"))
cg27193031_plot

## Mix and match for your own publication and based on your needs.






######################### Differentially Methylated Region analysis ###########################
## You might find in your research that there are genes implicated in whatever you are looking at
######################### Code for pulling out probes based on gene region ####################
## Search for Gene of interest on UCSC (need to have the 'EPIC.bed' file uploaded, i have this file!)
## I Have uploaded this file to https://github.com/pfransquet/EPIC-Bed-for-UCSC along with a short readme.
## Once on UCSC search for your gene of interest, i have put 7 genes of interest to one of my studies here for an example
## Zoom out a little to make sure you capture any nearby CpG Island
## Select a region (holding down control and sliding clicked mouse over) including some probes upstream and downstream, 
## if you highlight and select it will give you an exact measure of what you are looking at
## An example is for APP is around : chr21:27,152,684-27,736,037 - a huge 400kb region APP has lots of intronic regions
## in the code, name your subset data after the gene of interest, put in the chromosome e.g. ann.EPIC$chr=="chr9"
## then put in the location (i have minus and plus 1, just incase there is a probe just one bp outside the region)

## Sometimes you can see a number of probes on UCSC e.g. GAS1 gene shows 9 probes within the region, but when you put in the gene region for this code it might only 
## spit out 7, you can double check to make sure that your m file actually doesn't have data for the missing probes, in this case here is a code to see if there is actually data there
## probecheck <- subset(m, rownames(m)== "cg26447413")
## probecheck <- subset(m, rownames(m)== "cg02501714")
## For me it looks like these 2 probes cg26447413 and cg02501714 were removed in the probe removal process as part of the QC
## if you run code on probe that is there you will get a hit in your probecheck dataframe
## e.g probecheck <- subset(m, rownames(m)== "cg14416311")
## Because the probe is there it gives you data in the data frame

## Now to select all probes within a required gene region, and make a subset of the beta data frame.
## I have added gene identifiers to each of the probes so i know which gene they are in too, using 'rownames(beta_NR3C1) <- paste("NR3C1", rownames(beta_NR3C1), sep = "_")'

NR3C1_probes <- as.data.frame(subset(ann.EPIC, ann.EPIC$chr=="chr5" & ann.EPIC$pos> 142640231-1 & ann.EPIC$pos< 142851250+1))
NR3C1_probe_list <- NR3C1_probes$Name
beta_NR3C1 <- subset(beta, rownames(beta) %in% NR3C1_probe_list) 
rownames(beta_NR3C1) <- paste("NR3C1", rownames(beta_NR3C1), sep = "_")
rm(NR3C1_probes)
rm(NR3C1_probe_list)

NR3C2_probes <- as.data.frame(subset(ann.EPIC, ann.EPIC$chr=="chr4" & ann.EPIC$pos> 148912933-1 & ann.EPIC$pos< 149414545+1))
NR3C2_probe_list <- NR3C2_probes$Name
beta_NR3C2 <- subset(beta, rownames(beta) %in% NR3C2_probe_list)
rownames(beta_NR3C2) <- paste("NR3C2", rownames(beta_NR3C2), sep = "_")
rm(NR3C2_probes)
rm(NR3C2_probe_list)

FKBP5_probes <- as.data.frame(subset(ann.EPIC, ann.EPIC$chr=="chr6" & ann.EPIC$pos> 35501362-1 & ann.EPIC$pos< 35700206+1))
FKBP5_probe_list <- FKBP5_probes$Name
beta_FKBP5 <- subset(beta, rownames(beta) %in% FKBP5_probe_list)
rownames(beta_FKBP5) <- paste("FKBP5", rownames(beta_FKBP5), sep = "_")
rm(FKBP5_probes)
rm(FKBP5_probe_list)

BDNF_probes <- as.data.frame(subset(ann.EPIC, ann.EPIC$chr=="chr11" & ann.EPIC$pos> 27670032-1 & ann.EPIC$pos< 27747451+1))
BDNF_probe_list <- BDNF_probes$Name
beta_BDNF <- subset(beta, rownames(beta) %in% BDNF_probe_list)
rownames(beta_BDNF) <- paste("BDNF", rownames(beta_BDNF), sep = "_")
rm(BDNF_probes)
rm(BDNF_probe_list)

CRH_probes <- as.data.frame(subset(ann.EPIC, ann.EPIC$chr=="chr8" & ann.EPIC$pos> 67088056-1 & ann.EPIC$pos< 67091519+1))
CRH_probe_list <- CRH_probes$Name
beta_CRH <- subset(beta, rownames(beta) %in% CRH_probe_list)
rownames(beta_CRH) <- paste("CRH", rownames(beta_CRH), sep = "_")
rm(CRH_probes)
rm(CRH_probe_list)

CRHR1_probes <- as.data.frame(subset(ann.EPIC, ann.EPIC$chr=="chr17" & ann.EPIC$pos> 43853227-1 & ann.EPIC$pos< 43926705+1))
CRHR1_probe_list <- CRHR1_probes$Name
beta_CRHR1 <- subset(beta, rownames(beta) %in% CRHR1_probe_list)
rownames(beta_CRHR1) <- paste("CRHR1", rownames(beta_CRHR1), sep = "_")
rm(CRHR1_probes)
rm(CRHR1_probe_list)

CRHR2_probes <- as.data.frame(subset(ann.EPIC, ann.EPIC$chr=="chr7" & ann.EPIC$pos> 30682633-1 & ann.EPIC$pos< 30750929+1))
CRHR2_probe_list <- CRHR2_probes$Name
beta_CRHR2 <- subset(beta, rownames(beta) %in% CRHR2_probe_list)
rownames(beta_CRHR2) <- paste("CRHR2", rownames(beta_CRHR2), sep = "_")
rm(CRHR2_probes)
rm(CRHR2_probe_list)

## You can stick them all together into a beta file, as these genes are part of a cortisol analysis, i have called it cortisol_genes_beta
Cortisol_genes_beta <-rbind(beta_NR3C1,beta_NR3C2,beta_FKBP5,beta_BDNF,beta_CRH,beta_CRHR1,beta_CRHR2)
write.table(Cortisol_genes_beta, file="~/processing/Output/Cortisol gene probes beta.csv",sep=",",col.names=NA)


## I also want to make a couple of different data frames, one with both the targets file and the gene specific probe data
cort_genes <- t(Cortisol_genes_beta)
rownames(targets) <- targets$GID
targets_gene <-  merge(as.data.frame(targets), as.data.frame(cort_genes), by='row.names', all=TRUE)
View(targets_gene) ## the "Row.names" and your ID columns should line up.
targets_gene$Row.names = NULL
rownames(targets_gene) <- targets_gene$GID


## With my first analysis (a dementia data set) i exported and ran hundreds of ttests ANOVA and pwcorr in STATA, 
## Instead here is a way to to test distribution and run correlation as it is quicker to automate
install.packages("dplyr")
install.packages("Hmisc")
install.packages("ggpubr")
install.packages("psych")
install.packages("purrr")
install.packages("cowplot")

library("dplyr")
library("ggpubr")
library("Hmisc")
library("psych")
library("purrr")
library("cowplot")


## Can visualize single probe distribution with density plots e.g
ggdensity(targets_gene$NR3C1_cg17860381,main = "Density plot of NR3C1_cg21177852", xlab = "Methylation %")
## or a q-q plot which is sample vs normal distribution
ggqqplot(targets_gene$NR3C1_cg17860381)
## and can statistically check for nomality using the shapiro test (where significant results = non normal)
## This is also handy later on when looking at the linear models 
shapiro.test(targets_gene$NR3C1_cg17860381)


## but what we really want to do is see if all probes are normal, so we need to automate it!
## I found this function here "https://stackoverflow.com/questions/33489330/how-to-test-the-normality-of-many-variables-in-r-at-the-same-time"
## First just have the probes of interest from your targets_gene file
Probes_only <- as.data.frame(targets_gene[,c(40:510)]) 
## Important here that this file is a data frame or the process wont work, also best to use 'targets_gene' which has the cortisol and the ID matched up
## this will keep them in check for the next bit

cortprobelist <- as.character(colnames(Probes_only))

shapiro_test_df <- function(df, bonf= TRUE, alpha= 0.05) {
  l <- lapply(df, shapiro.test)
  s <- do.call("c", lapply(l, "[[", 1))
  p <- do.call("c", lapply(l, "[[", 2))
  if (bonf == TRUE) {
    sig <- ifelse(p > alpha / length(l), "H0", "Ha")
  } else {
    sig <- ifelse(p > alpha, "H0", "Ha")
  }
  return(list(statistic= s,
              p.value= p,
              significance= sig,
              method= ifelse(bonf == TRUE, "Shapiro-Wilks test with Bonferroni Correction",
                             "Shapiro-Wilks test without Bonferroni Correction")))
}

options(scipen=4)
options(digits=4)
testresults <-as.data.frame (shapiro_test_df(Probes_only)) ## H0 will be normal where Ha will be non-normal
testresults$probe <- cortprobelist
rownames(testresults) <- NULL

## Now make seperate data frames for normal and non-normally distributed probes

nonnormalrow <- as.character(which(testresults$p.value <= 0.05))
nonnormal <- as.character(testresults[nonnormalrow,"probe"])

normalrow <- as.character(which(testresults$p.value > 0.05))
normal <- as.character(testresults[normalrow,"probe"])

casd <- as.character (paste("Cortisol"))
nonnormalcort <- c(casd,nonnormal)
normalcort <- c(casd,normal)

## Now that we have lists of which probes are normally and non normally distributed, we can run parametric (pearsons) and non parametric (spearmans)
## Make a dataframe with the probe data and cortisol measures 
cormat <- targets_gene[,c(30,40:510)] ## Make sure the cortisol measure matches the rowname ID

## now subset based on normality of data (made from previous lists)
nonnormal_probes <- cormat[,nonnormalcort ]
normal_probes <- cormat[,normalcort ]

## Correlate the non normal probes using spearmans method, adjusting for multiple testing using holm
corrtestholm <- corr.test(nonnormal_probes, method="spearman", adjust="holm")

r <- as.data.frame(corrtestholm[["r"]])
p <- as.data.frame(corrtestholm[["p"]])
se <- as.data.frame(corrtestholm[["se"]])

names(r)[names(r) == "Cortisol"] <- "r"
r$p.val.holm <- p$Cortisol
r$se <- se$Cortisol

## Now make a results table with the correlation to continuous measure, the adjusted pval and the standard error
ncol(r)
## Your correlation to the continuous measure should be in the first column, and
## your p val and se will be the last columns in the 'r' data table, for me it is 188, 189.  
nonnorm_corresults <- r[,c(1,188,189)]

## Now check correlation for the normal probes
corrtestholm <- corr.test(normal_probes, method="pearson", adjust="holm")

r <- as.data.frame(corrtestholm[["r"]])
p <- as.data.frame(corrtestholm[["p"]])
se <- as.data.frame(corrtestholm[["se"]])

names(r)[names(r) == "Cortisol"] <- "r"
r$p.val.holm <- p$Cortisol
r$se <- se$Cortisol

ncol(r) ## remember pick the first (r) and last two columns (p.val and se)

norm_corresults <- r[,c(1,287,288)]


## Can stick the correlation tables together! remove the top row, dont need it anyway, just shows how cortisol correlates with its self
nonnorm <- nonnorm_corresults [-c(1),]
norm <- norm_corresults[-c(1),]
All_Corr_Results <- rbind(nonnorm,norm)
All_Corr_Results$probe <- rownames(All_Corr_Results)

## the amount of observations should match that of the variables from the "Probes_only" e.g. there are data for all 471 probes that i pulled out from the gene regions.
## Now pull out the data for those that passed p.val.adj. (i used the above example of making a list of rows, then pulling out the probes,
## but it didnt like the "as.character". i dont know why, but when removing it and using the list as an integer it worked

corrprobe_rowlist <- which(All_Corr_Results$p.val.holm <0.05)
sigcorrprobes <- All_Corr_Results[corrprobe_rowlist,]
write.csv(sigcorrprobes, file="~/processing/Output/sigcorrprobes.csv")


## 23 probes had methylation correlated with cortisol levels. Now lets graph it. I want to try to automate it with a function as ggplot requires graphs to be made one by one,
## but we want the same graph with different data so its best if we try and automate
sigcorrprobe_list <- as.character(sigcorrprobes$probe)
corrgraphtable <- Probes_only[,sigcorrprobe_list]
corrgraphtable$Cortisol <- targets_gene$Cortisol ## Just check along the way that data matches with the ID, as sometimes when merging data sets and manipulating data, the order changes. so data doesnt line up


## The question is whether methylation explains difference in cortisol, therefor methylation is the explanatory variables and cortisiol the response variable
## https://aosmith.rbind.io/2018/08/20/automating-exploratory-plots/
expl = names(corrgraphtable)[1:23] ## the sig probes
response = names(corrgraphtable)[24] ## Continuous measure 

response = set_names(response)
response

expl = set_names(expl)
expl

## EXAMPLE FUNCTION
#scatter_fun = function(x, y) {
#  ggplot(dat, aes(x = .data[[x]], y = .data[[y]]) ) +
#    geom_point() +
#    geom_smooth(method = "loess", se = FALSE, color = "grey74") +
#    theme_bw()
#}

## My modified function for this code
scatter_fun = function(x, y) {
  ggplot(corrgraphtable, aes(x = .data[[x]], y = .data[[y]]) ) +
    geom_point(size=2) +
    geom_smooth(method =lm, se = FALSE, color = "red") +
    theme_bw()
}


corr_plots = map(expl, ~scatter_fun(.x, "Cortisol") )


pdf("Correlation Plots_All Data_04OCT21.pdf")
corr_plots
dev.off()


## I want to see the correlation when compared between a binary variable (here for example maternal pregnancy PTSD groups), or catagorical variable (Low, normal and high cortisol) 
## you can add it into the corrgraphtable and then plot them both
corrgraphtable$Preg_PTSD <- targets_gene$Pregnancy_PTSD
corrgraphtable$CortCat <- targets_gene$TertileCort
## If numbers ggplot will think that the colour variable is continuous so renaming helps
corrgraphtable$Preg_PTSD_Cat <- ifelse(corrgraphtable$Preg_PTSD=="0", "No", ifelse(corrgraphtable$Preg_PTSD=="1", "Yes", NA))  
corrgraphtable$CortCat_cat <- ifelse(corrgraphtable$CortCat=="0", "Normal", ifelse(corrgraphtable$CortCat=="1", "Low", ifelse(corrgraphtable$CortCat=="2", "High", NA)))  


scatter_fun = function(x, y) {
  ggplot(corrgraphtable, aes(x = .data[[x]], y = .data[[y]], color = Preg_PTSD_Cat) ) +
    geom_point(size=2) +
    geom_smooth(method =lm, se = FALSE, fullrange=TRUE) +
    theme_bw()
}


corr_plots = map(expl, ~scatter_fun(.x, "Cortisol") )


pdf("Correlation plots, PTSD Seperarted_04OCT21.pdf")
corr_plots
dev.off()


scatter_fun = function(x, y) {
  ggplot(corrgraphtable, aes(x = .data[[x]], y = .data[[y]], color = CortCat_cat ) ) +
    geom_point(size=2) +
    geom_smooth(method =lm, se = FALSE, fullrange=TRUE) +
    theme_bw()
}


corr_plots = map(expl, ~scatter_fun(.x, "Cortisol") )


pdf("Correlation plots, Cortisol Catagory Seperarted_04OCT21.pdf")
corr_plots
dev.off()



######################### REGRESSION OF SIGNIFICANT PROBES IN CANDIDATE GENE ANALYSIS ############################
## Make a table with m values (as this is better for analysis) from your data of interest, you can use sig list from correlation analysis.
## Beware that some of this code might over write some of the limma model code from above due to naming conventions
sigcorrprobes <- read.csv(file="~/processing/Output/sigcorrprobes.csv",row.names = 1)
library(stringr)
library(limma)
probesplit <- as.data.frame (str_split_fixed(sigcorrprobes$probe, "_", 2))
sigcorrprobes$probes <- probesplit$V2
View(sigcorrprobes)

sigprobes <- as.character(sigcorrprobes$probes)
sigcorrprobes$probes

m_sig <- as.data.frame(m[sigprobes,])
beta_sig <- as.data.frame(beta[sigprobes,])

## Lets work on some multivariate linear regression
## First want to check the distribution of methylation in the sig probes
## This uses some data already made above 
## This scatter plot will take the significant probes from the continuous cortisol and compare methylation between PTSD vs Non-PTSD
pdf(file="Sig 23 Scatter Box Plots.pdf", width = 14)
par(mfrow=c(2,5))
#write.csv(bl_list, file = "A:/Data/R Analysis No SNP/DMPs/BL Top 200 dec.csv")
plotCpg(dat=m,pheno=targets_gene$Cortisol, cpg=sigprobes, type = "continuous", measure = "beta",
        ylab = "M", fitLine = TRUE) ## plotCpg is a function of minfi
dev.off()

## Now lets do some linear modelling 
install.packages("tidyverse")
library(tidyverse)
library(ggpubr)

## Lets do a test model with the first sig CPG
sigprobes
## for me it is  cg17860381, which in my siggcorrporbes table i can see is an NR3C1 probe
testmodel <- lm(Cortisol ~ NR3C1_cg17860381 , data = targets_gene) 
summary(testmodel)
## From previous corr test r should be -0.2242, and here i got r squared of 0.0282 (r = -0.167).
## Might be a bit different as the methylation is not normally distributed.
## Remember you can visualise and even stat check normality 
## Density plot
ggdensity(targets_gene$NR3C1_cg17860381,main = "Density plot of NR3C1_cg21177852", xlab = "Methylation %")
## q-q plot (sample vs normal distribution)
ggqqplot(targets_gene$NR3C1_cg17860381)
## and can statistically check for normality using the Shapiro test (where significant results = non normal)
shapiro.test(targets_gene$NR3C1_cg17860381)

## This one may need a nonparametric regression.
## I would consult a statistician before going down this path but here is a package you can use if wanted
## http://users.stat.umn.edu/~helwig/notes/smooth-notes.html
shapiro.test(targets_gene$CRHR2_cg09797340)



## Lets do a regression on a normally distributed probe, we will go with CRHR2_cg09797340
## First check normality
shapiro.test(targets_gene$CRHR2_cg09797340)
## Now linear model
testmodel <- lm(Cortisol ~ CRHR2_cg09797340 , data = targets_gene) 
summary(testmodel)
## From previous corr test r should be -0.2307, and here i got 0.0532 (r squared), which is spot on (r = -0.23065)
## Can add other variables to the model, lets use what we had in the Limma model, Age_Child, Gender_Child, Age_Mother, Chip, counts.CD8T,
## counts.CD4T, counts.NK, counts.Bcell, counts.Mono and counts.Neu.
## This time rather than remove the samples that have maternal age missing, i will quickly impute them using mice
## Also will only impute numbers so if you need to impute categories which are text, will have to convert to numbers first
install.packages("mice")
library(mice)
targ <- mice(targets_gene, m = 1, pri = FALSE, maxit = 0, seed = 75631564)
targ
targets_gene_imp <- complete(targ,1)
View(targets_gene_imp)

## Also lm for modelling doesnt like '.' in variable names

## now the model
testmodel <- lm(Cortisol ~ CRHR2_cg09797340 + Age_Child + as.factor(Gender_Child) + Age_Mother  + as.factor(Chip)  +
                counts.CD8T  +  counts.CD4T  + counts.NK  + counts.Bcell  + counts.Mono  + counts.Neu, data = targets_gene_imp) 
summary(testmodel)

## With this model the relationship between methylation and cortisol is lost







#################################################################################################################################

######################################################### DMRcate  ##############################################################
library(DMRcate)
## You can use "DMRcate" to look for differentially methylation regions (DMRs) across the genome.
## It used to work regardless but i think since an update it only reports if it has findings, so here is the code i used when it did work but cant check it as i have no 
## significant DMRs in my data set.

## First Need to set up a design matrix just like limma
## This is a rebuild of covar4 from above, but i only need one line to name, because its already built
Cortisol <- targets_gene$Cortisol
designcovar <- model.matrix(~Cortisol + targets_gene$Age_Child + as.factor(targets_gene$Gender_Child) + as.factor(targets_gene$Chip))
head(designcovar)

## If you dont want to include CpG's that werent significant between groups, i wouldnt suggest this as DMR cate gives weights to 
## probes and you might miss a good region just because some non differential CpGs are supposed to be present.
## you can  make a subset of 'm' and 'beta' based on the top limma probes (non adjusted pval <.05)

## Use cpg.annotate to let the code know where probes are situated on each gene
ann.DMR_CORT <- cpg.annotate("array", m, what="M", arraytype = "EPIC", 
                             analysis.type="differential", design=designcovar, fdr = 0.05, coef=2)

## If nothing is found you will get the error:
## Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.
## So for this example i will raise the FDR so we have something to work with 
## remember you can use for example '??cpg.annotate' which if the library is installed should generatlly 

ann.DMR_CORT <- cpg.annotate("array", m, what="M", arraytype = "EPIC", 
                             analysis.type="differential", design=designcovar, fdr = 0.1, coef=2)

## for me, there was still nothing at fdr = 0.1, so we can use the following code to look at the pval range and set a pval different to include more CpG's (but increases type I error)
DMRpval <-as.matrix(ann.DMR_CORT@ranges@elementMetadata@listData[["ind.fdr"]])
DMRpvaldf <- as.data.frame(DMRpval)
summary(DMRpvaldf)

## So basically in my data set min FDR is 0.410 so lets use 0.5 and see what happens (not recommended, just for example)
ann.DMR_CORT <- cpg.annotate("array", m, what="M", arraytype = "EPIC", 
                             analysis.type="differential", design=designcovar, fdr = 0.5, coef=2)

## This time i get the notice 
## "Your contrast returned 29684 individually significant probes. We recommend the default setting of pcutoff in dmrcate()."

## Now make the DMR file, you can change the p value here too which will allow more DMR but increases type I error
## Lambda refers to how many base pairs max between each significant probe
DMR_CORT_output <- dmrcate(ann.DMR_CORT, lambda=1000, pcutoff = 0.05, C=2)
## You can also save this, and name as per model you used to get it
save(DMR_CORT_output,file="~/processing/R Data/DMRs/DMR_PTSD_output_age_sex_chip_04OCT2021.Rdata")

## I think you might need the internet for this next bit
results.ranges_DMR_CORT_output <- extractRanges(DMR_CORT_output, genome = "hg19")

## DMRcate doesnt allocate a DMR number to the DMR, it just goes down the results ranges list when you choose which DMR to plot.
## If you want to go in and look at the data and give the DMR a number do the following

RR_CORT <- as.data.frame (results.ranges_DMR_CORT_output)
View(RR_CORT)
RR_CORT$DMRID <- row.names(RR_CORT)
row.names(RR_CORT) =NULL
RR_CORT$DMRNO <- row.names(RR_CORT)
row.names(RR_CORT) = RR_CORT$DMRID
RR_CORT$DMRID = NULL
## You can order by whichever value you want, i will order here by Fisher’s multiple comparison statistic
RR_CORT <- RR_CORT[order(RR_CORT$Fisher), , drop = FALSE]
View(RR_CORT)

## Now you can look into the RR dataframe, and plot a particular DMR.
## First we will subset the results so we can prioritise Fisher’s multiple comparison statistic, and then also the number of CpGs we want in a DMR.
summary(RR_CORT$Fisher)
## Lets use a very non-conservative multiple testing cutoff
RR_CORT.top <- subset(RR_CORT,Fisher< 0.7)

## Now lets set a minimum CpG amount. e.g DMR with at least 3 probes
RR_CORT.top <- subset(RR_CORT.top,no.cpgs>=3)
## Now we have a smaller table of more interesting results.
## You can use other cut offs too if you want 
View(RR_CORT.top)
write.table(RR_CORT.top, file="~/processing/R Data/DMRs/DMR_CORT_top_DMR.csv",sep=",",col.names=NA)

## To list what probes are in which DMR (this is really crude but it works!)
cgID_CORT <- as.data.frame(dmrcoutput.adj.bl$input)
#cgID$keep <- as.numeric(cgID$step.dmr=="TRUE")
#cgID <- subset(cgID,keep=="1")

#To get the beta's from the probes in a particular DMR, choose a DMR from the list and pull out the betas like we did the single probes
DMR_629_probes <- as.data.frame(subset(ann.EPIC, ann.EPIC$chr=="chr1" & ann.EPIC$pos> 153362927-1 & ann.EPIC$pos< 153364020+1))
DMR_629_probe_list <- DMR_629_probes$Name
beta_DMR_629 <- subset(beta, rownames(beta) %in% DMR_629_probe_list) 
rownames(beta_DMR_629) <- paste("DMR_629", rownames(beta_DMR_629), sep = "_")
View (beta_DMR_629)
rm(DMR_629_probes)
write.table(beta_DMR_629, file="~/processing/R Data/DMRs/DMR_629.csv",sep=",",col.names=NA)

## Can make a scatter plot
pdf(file="~/processing/R Data/DMRs/DMR_629.pdf", width = 14)
par(mfrow=c(1,4))
plotCpg(dat=beta,pheno=targets_gene$Cortisol, type="continuous", cpg=DMR_629_probe_list, ylab = "Beta Values")
dev.off()

## you can also do categorical plots
pdf(file="~/processing/R Data/DMRs/DMR_629_cat.pdf", width = 14)
par(mfrow=c(1,4))
plotCpg(dat=beta,pheno=targets_gene$TertileCort, type="categorical", cpg=DMR_629_probe_list, ylab = "Beta Values")
dev.off()

## DMRcate also has a plotting code within it, but it can get quite messy if you have more than n=12
## Also this only works if connected to the internet
group <- c('0'="magenta", '1'="royal blue")
cols <- group[as.character(targets_gene$Pregnancy_PTSD)]

DMR.plot(ranges=results.ranges_DMR_CORT_output, dmr=2, CpGs=beta, what="Beta", arraytype = "EPIC", phen.col=cols, genome="hg19")



################################################### Gene ontology using 'gometh'#########################################################
## There are other ways to do pathway analysis but this is the one i used.
## Go meth is a missMethyl function
## library (missMethyl)
## With gometh you need a set of significant results, and then a set of all probes to compare against
## Lets take the non-multiple testing significant results from the first limma model
sig.cpg <- subset(limmatable_covarRun1, P.Value <= 0.05)
sig.cpg$probe <- rownames(sig.cpg)
sig.cpg <- as.character(sig.cpg$probe)

## Now take the whole list of CpGs, from the gmset with all the problem probes removed from the very start of the pipeline
all.cpg <- as.character(gmset.pr@rowRanges@ranges@NAMES)

## you can save these lists if you want too
save(sig.cpg,file="~/processing/R Data/GOMETH/sig.cpg.Rdata")
save(all.cpg,file="~/processing/R Data/GOMETH/all.cpg.Rdata")

## Or make a list to look at in excel
write.csv(sig.cpg,file="~/processing/R Data/GOMETH/sig.cpg.csv")
write.csv(all.cpg,file="~/processing/R Data/GOMETH/all.cpg.csv")

################################################## Next part needs to be done with internet access #######################################
## There are two different pathway data bases you can look at, KEGG and GO. 
## First KEGG
par(mfrow=c(1,1))
KEGG <- gometh(sig.cpg, all.cpg = all.cpg, collection = c("KEGG"), array.type = c("EPIC"), plot.bias = TRUE, prior.prob = TRUE)
save(KEGG,file="~/processing/R Data/GOMETH/KEGG_PTSD.Rdata")
write.table(KEGG, file="~/processing/KEGG_PTSD.csv",sep=",",col.names=NA)


## And GO
par(mfrow=c(1,1))
GO <- gometh(sig.cpg, all.cpg = all.cpg, collection = c("GO"), array.type = c("EPIC"), plot.bias = TRUE, prior.prob = TRUE)
save(GO,file="~/processing/R Data/GOMETH/GO_CORT.Rdata")
write.table(GO, file="~/processing/GO_CORT.csv",sep=",",col.names=NA)



