# Genomewide association analysis (GWAS)

## Objectives
In this practical, you will learn the basic steps of data quality assessment and QC on genotype data that are typically carried out in genome-wide association studies (GWAS).
* [Genotype Data format](#genotype-data-format)
* [Sample QC](#Sample-QC)
* [Variant QC](#Variant-QC)
* [Association analysis](#Assoc)

## Prerequisites
For sample and variant QC, we need the following software that is pre-installed in the VM
- `PLINK`[1.9 beta](https://www.cog-genomics.org/plink/) for data management 
- `R` for visualization of different QC measures
  
## Step 0: Download the genotype data


- Download all zip file under it to your current directory
```bash
wget https://github.com/claratsm/HTI5802_Algorithms-in-Bioinformatics-and-Genomics/Practical/gwas.tar.gz
tar -zxvf gwas.tar.gz
```

- Move the files in 1000gData to your current directory !!!!!!!!!!!!!!!!!!!!!!!! NEED TO MAKE THE DATA SMALLER FIRST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```bash
mv ~/1000gData/* .
```

- First, go to a directory named `gwas` for this practical section
```bash
cd ~/gwas
```

Use `ls` to list the files in your current directory. You should get `NINE` files, including the two sets of PLINK binary genotype files (`BED`, `BIM` and `FAM`)
- Now your current folder should contain these files
> ASA.1000G.to-update-name.snp<br>
> chrAll.ASA.1000GP-All.bed<br>
> chrAll.ASA.1000GP-All.bim<br>
> chrAll.ASA.1000GP-All.fam<br>
> chrAll.ASA.bed<br>
> chrAll.ASA.bim<br>
> chrAll.ASA.fam<br>
> integrated_call_samples_v3.20130502.ALL.panel<br>
> practical2.PCAplot.R<br>

### Dataset
The dataset used in this practical was simuated from haplotypes of East Asian samples of the 1000 Genomes Project ([Phase 3](https://www.internationalgenome.org/category/phase-3/)). SNPs included in the dataset reflect to those assayed in the [Illumina Asian Screening array](https://www.illumina.com/products/by-type/microarray-kits/infinium-asian-screening.html) designed to maximize the genomic coverage for East Asian population.

## Genotype Data format
### - PLINK binary variant file (BIM)
Let's have a look at the variant `BIM` file
```bash
head chrAll.ASA.bim
```
BIM file stores the variant information and must contain as many markers as are in the BED file.<br>
It contains no header line and each line describes a single marker with 4 columns:

> 1. Chromosome code (1-22, 23/X, 24/Y, 25/XY, 26/MT; '0' indicates unknown)<br>
> 2. Variant identifier<br>
> 3. Position in morgans or centimorgans (safe to use dummy value of '0')<br>
> 4. Base-pair coordinate (1-based)
> 5. A1 allele (Default: minor allele)
> 6. A2 allele (Default: major allele)

:green_book: **Q:** How many SNPs are there?
```bash
wc -l chrAll.ASA.bim
```
:green_book: **Q:** Where are these SNPs located in?
```bash
cut -f 1 chrAll.ASA.bim | sort | uniq -c
```

### - PLINK binary pedigree file (FAM)
Let's have a look at the `FAM` file
```bash
less -S chrAll.ASA.fam   # type 'q' to quit
```
The FAM file is a white-space (space or tab) delimited file<br>
It has no header line, and one line per sample with 6 columns:

> 1. Family ID ('FID')<br>
> 2. Individual ID('IID')<br>
> 3. Within-family ID of father ('0' if father isn't in dataset)<br>
> 4. Within-family ID of mother ('0' if mother isn't in dataset)<br>
> 5. Sex (1=male; 2=female; 0=unknown)<br>
> 6. Phenotype (1=control, 2=case, -9 / 0=missing for case-control; numeric for quantitative trait)<br>

:green_book: **Q:** How many samples are there?
```bash
wc chrAll.ASA.fam
```
:green_book: **Q:** Any sample belongs to the same family?
```bash
awk '$3!=0 || $4!=0' chrAll.ASA.fam
```

## Sample-QC

## Why do we need the QC on GWAS data?
Suboptimal study design, poor sample quality, and errors in genotype calling can introduce systematic biases into GWAS, leading to spurious associations. A thorough QC can help us identify samples and variants that should be removed prior to association analysis to minimize false-positive and false-negative associations.

Here we assume that the study design has been conducted appropriately and basic QC to genotypes (e.g. QC based on GenTrain score in Illumina GenomeStudio or rare SNPs calling using zCall) has been performed after clustering and calling from probe intensity data. 

## The QC protocol: from Sample QC to Variant QC 
The QC protocol for GWAS typically split into two broad categories, `Sample QC` and `Variant QC`. Sample QC is usually done prior to Variant QC in order to maximise the number of variants remained for association analysis.

<p align="center">
<img src="https://user-images.githubusercontent.com/8644480/170681398-e29f945e-fc94-4876-b695-9a8f2250968e.png"  width="500" height="300">
</p>

In this practical, we will investigate the
+ **Missingness** for samples and variants
+ **Heterogeneity** in terms of inbreeding coefficient for sex chromosomes and autosomes
+ Identify-by-descent (IBD) and coefficient of relationship (PI_hat) for **relatedness** assessment
+ Principal component analysis (PCA) for identification of **population outliers**
+ Hardy-weinberg equilibrium
+ 
## Step 1: Individuals with excessive missing genotypes
- Obtain profile of missingness per individual and per SNP
```bash
~/plink --bfile chrAll.ASA --missing --out chrAll.ASA.beforeQC
```
This command generates two files, including `.imiss` for sample-based and `.lmiss` for variant-based missingness
> chrAll.ASA.beforeQC.imiss    # per individual<br>
> chrAll.ASA.beforeQC.lmiss    # per variant<br>

For both files, the last three columns measure the missingness for each individual or variant 
> N_MISS: Number of missing genotype call(s)<br>
> N_GENO: Number of genotype call(s)<br>
> F_MISS: Missing call rate<br>

#### **QC: Record the sample IDs with call rate <= 0.98 for removal**
```bash
awk 'NR>1 && $6>0.02 { print $1,$2 }' chrAll.ASA.beforeQC.imiss > to-remove.mind02.indiv
cut -d' ' -f 2 to-remove.mind02.indiv > to-remove.mind02.iid
```

- Plot the distribution of missingness 
We can then use R script to generate histogram of missingness
```R
# ========================== R code ==========================
# ----------------------------------------------------------#
#         (1)  SAMPLE CALL RATE    - threshold = 98%        #
# ----------------------------------------------------------#
imiss <- read.table("chrAll.ASA.beforeQC.imiss",h=T)
head(imiss)
summary(imiss$F_MISS)

# Plot missingness across samples
hist(imiss$F_MISS, freq=T, col="darkred", border ="black", main="Sample Call Rate", 
xlab="F_MISS", ylab="Number of samples")
```
```R
# Plot missingness with altered y-axis for a zoom in view
hist(imiss$F_MISS, breaks=seq(0,0.6,0.05), freq=T, col="darkred", border ="black", 
main="Sample Call Rate", xlab="F_MISS", ylab="Number of samples",ylim=c(0,10))
# ============================================================
```
![practical2 missing-hist1](https://user-images.githubusercontent.com/8644480/170830897-4eb7a450-bc51-42ea-b68b-bf5ea98197ab.png)
![practical2 missing-hist2](https://user-images.githubusercontent.com/8644480/170830900-ef75c54a-ba80-4393-9c77-7502ac7c319a.png)

:closed_book: **Q:** Can you try to plot the right number of excluded samples?
<details>
  <summary> Try your own R codes </summary>
<p></p>

- **Answer 1**
```R
#==== R =====
hist(imiss$F_MISS, breaks=seq(0,0.6,0.01), freq=T, col="darkred", border="black", main="Sample Call Rate", xlab="F_MISS", ylab="Number of samples", ylim=c(0,100))
abline(v=0.02, lwd=2, lty=2, col="darkblue")
abline(v=0.01, lwd=2, lty=2, col="darkgreen")
```
![practical2 missing-hist3](https://user-images.githubusercontent.com/8644480/170831053-f941b2df-090a-4880-8711-0dd7a4b3fb77.png)

- **Answer 2**
```R
#==== R =====
plot(sort(imiss$F_MISS), pch=20, col="darkred", main="Sample Call Rate", xlab="ASA samples", ylab="F_MISS")
abline(h=0.02, lwd=2, lty=2, col="darkblue")
abline(h=0.01, lwd=2, lty=2, col="darkgreen")
```
![practical2 missing-hist4](https://user-images.githubusercontent.com/8644480/170831061-7d3811a0-8502-4a6e-bc22-e4d5a45eb52e.png)

- **Answer 3**
```R
#==== R =====
plot(imiss$F_MISS,pch=20,col="darkred", main="Sample Call Rate", xlab="ASA samples", ylab="F_MISS")
abline(h=0.02, lwd=2, lty=2, col="darkblue")
abline(h=0.01, lwd=2, lty=2, col="darkgreen")
```
![practical2 missing-hist5](https://user-images.githubusercontent.com/8644480/170831069-58dd61a8-b4db-4ee2-a6a7-e63968b35108.png)

</details>
We typically remove samples and variants with high degree of missingness using less stringent threshold (e.g. missingness > 5% or 10%) at early stages of the QC steps; however, to illustrate how bad calling related to different QC measures, we keep all samples till the end of the QC practical.

## Step 2: Individuals with sex discrepancy
- Obtain missingness of chr X
```bash
~/plink --bfile chrAll.ASA --chr 23 --missing --out chrX.ASA.beforeQC
```
- Check sex
```bash
~/plink --bfile chrAll.ASA --check-sex --out chrX.ASA
```
The function of `--check-sex` normally compares the sex assignments in the input pedigree file (`PEDSEX`) with inbreeding coefficients (F) imputed from SNPs on chromosome X. By default, the F estimates smaller than 0.2 yield female calls, and values larger than 0.8 yield male calls for `SNPSEX`. 

```R
# ========================== R code ==========================
sexcheck <- read.table("chrX.ASA.sexcheck",h=T)
head(sexcheck)

#       FID      IID PEDSEX SNPSEX  STATUS       F
# 1 id1_0211 id2_0211      2      2      OK -1.0180
# 2 id1_1202 id2_1202      2      0 PROBLEM  0.6545
# 3 id1_1130 id2_1130      2      2      OK -0.7449
# 4 id1_1064 id2_1064      2      2      OK -0.2048
# 5 id1_1109 id2_1109      2      0 PROBLEM  0.7954
# 6 id1_0108 id2_0108      2      0 PROBLEM  0.3521

mismatch <- sexcheck[sexcheck$STATUS=="PROBLEM",]

colsex <- c("darkblue","darkred")
plot(sexcheck$F, col=colsex[sexcheck$PEDSEX], main="Sex check", pch=20, xlab="ASA samples", ylab="chrX Inbreeding coefficient (F)")
points(row.names(mismatch),mismatch$F,pch=22,bg="green",col="black",lwd=2,cex=2)
abline(h=0.2, lwd=2, lty=2, col="red")
abline(h=0.8, lwd=2, lty=2, col="blue")
legend("bottomright",c("Male PEDSEX","Female PEDSEX","sample with PROBELM"), col=c(colsex,"black"),pt.bg="green", pch=c(20,20,22))
# =============================================================
```
![practical2 sexcheck](https://user-images.githubusercontent.com/8644480/170831468-747d54ae-1b5a-4057-aa00-319f9a2a4864.png)

- What is the relationship between missingness and inbreeeding coefficient for chrX?
```R
# ========================== R code ==========================
imiss.X <- read.table("chrX.ASA.beforeQC.imiss",h=T)
sexcheck.imiss <- merge(imiss.X, sexcheck, by="IID")
mismatch.imiss <- sexcheck.imiss[sexcheck.imiss$STATUS=="PROBLEM",]
plot(sexcheck.imiss$F_MISS, sexcheck.imiss$F, pch=20, col=colsex[sexcheck.imiss$PEDSEX], xlab="Number of missing genotypes on chrX", ylab="chrX Inbreeding coefficient (F)")
points(mismatch.imiss$F_MISS, mismatch.imiss$F, pch=22, bg="green", col="black", lwd=2, cex=1.5)
abline(h=0.2, lwd=2, lty=2, col="red")
abline(h=0.8, lwd=2, lty=2, col="blue")
legend("bottomleft", c("Male PEDSEX","Female PEDSEX","sample with PROBELM"), col=c(colsex,"black"),pt.bg="green", pch=c(20,20,22))
# =============================================================
```
![practical2 sexcheck-fmiss](https://user-images.githubusercontent.com/8644480/170832150-fb1471ed-be27-4e1e-80eb-cdf78989d21f.png)

In addition to poor sample or genotyping quality, X chromosome aneuploidy, such as Turner syndrome (e.g. 45,X0) and Klinefelter syndrome (47, XXY), may lead to abnormal heterogenity. Missingness for SNPs on chr Y can be used to impute sex using `--check-sex ycount [female max F] [male min F] [female max Y obs] [male min Y obs]`. Before determining the minimum and maximum threshold of observed chrY variants, we can first set [female max Y obs] to maximum number of chrY variants and [male min Y obs] to 0  
```bash
~/plink --bfile chrAll.ASA --check-sex ycount 0.2 0.8 805 0 --out chrXY.ASA.without-ycount
```
:closed_book: **Q:** Do the samples detected by sex discrepancy or abnormal heterogenity look like having X chromosome aneuploidy?
<details>
  <summary> Try your own PLINK / R codes </summary>
<p></p>
  
- ** 1 **
```bash
egrep PROBLEM chrXY.ASA.without-ycount.sexcheck
#         FID       IID PEDSEX SNPSEX  STATUS       F YCOUNT
#    id1_1202  id2_1202      2      0 PROBLEM 0.65450    607
#    id1_1109  id2_1109      2      0 PROBLEM 0.79540    636
#    id1_0108  id2_0108      2      0 PROBLEM 0.35210    569
#   CHSHet002 CHSHet002      1      0 PROBLEM 0.67790    804
#     id1_300   id2_300      1      2 PROBLEM 0.03826      2
#     id1_301   id2_301      2      1 PROBLEM 1.00000    799
#     id1_500   id2_500      1      2 PROBLEM 0.01425      2 
#     id1_501   id2_501      2      1 PROBLEM 1.00000    805
```
- ** 2 **
```R
# ========================== R code ==========================
sexcheck.XY <- read.table("chrXY.ASA.without-ycount.sexcheck",h=T)
mismatch.XY <- sexcheck.XY[sexcheck.XY$STATUS=="PROBLEM",]
colsex <- c("darkblue","darkred")
plot(sexcheck.XY$YCOUNT, sexcheck.XY$F, pch=20, col=colsex[sexcheck.XY$PEDSEX], xlab="Number of non-missing genotypes on chrY", ylab="chrX Inbreeding coefficient (F)")
points(mismatch.XY$YCOUNT, mismatch.XY$F, pch=22, bg=colsex[mismatch.XY$PEDSEX], col="white", lwd=1, cex=1.5)
abline(h=0.2, lwd=2, lty=2, col="red")
abline(h=0.8, lwd=2, lty=2, col="blue")
# =============================================================
```  
![practical2 sexcheck-XY](https://user-images.githubusercontent.com/8644480/170838631-9463a53c-2f21-4151-bf44-b91d5d33111e.png)
</details> 
  
You can now rerun `--check-sex` with appropriate [female max Y obs] and [male min Y obs] thresholds, e.g. 100 and 700. Samples with abormal number of non-missing genotypes on chrY will be flagged as ambiguous in `SNPSEX`.         
```bash
~/plink --bfile chrAll.ASA --check-sex ycount 0.2 0.8 100 700 --out chrXY.ASA.with-ycount
egrep -h PROBLEM chrXY.ASA.without-ycount.sexcheck chrXY.ASA.with-ycount.sexcheck | sort | uniq -u
```
:closed_book: **Q:** What is the most likely karyotype for the **id2_669** sample?
<details>
  <summary> Answer </summary>
  
 The **id2_669** individual is likely to have Turner syndrome with 45,X0 or 45,X karyotype. Copy number variation analysis on SNP array data is needed to confirm the number of copy of X chromosome.
</details>
  
  
#### **QC: Record the sample IDs with sex discrepancy problem for removal**
```bash
awk '$5=="PROBLEM"{ print $1,$2 }' chrXY.ASA.with-ycount.sexcheck > to-remove.sexmismatch.indiv
```

## Step 3: Individuals with outlying heterozygosity rate
To avoid bias by genotyping error of rare variants and linkage disequilibrium, we usually perform the heterogeneity check using only common variants (MAF>=5%) and SNPs in strong LD
- Filter out low frequency and rare variants with MAF cut-off of 5%
```bash
~/plink --bfile chrAll.ASA --autosome --maf 0.05 --make-bed --missing --out chr1-22.ASA.maf05
```
- LD pruning using R-squared 0.1
```bash
~/plink --bfile chr1-22.ASA.maf05 --indep-pairwise 200 50 0.1 --out chr1-22.ASA.maf05.pruning
```
```bash
~/plink --bfile chr1-22.ASA.maf05 --extract chr1-22.ASA.maf05.pruning.prune.in --make-bed --out chr1-22.ASA.maf05.pruned
```
```bash
~/plink --bfile chr1-22.ASA.maf05.pruned --het --out chr1-22.ASA.maf05.pruned
```
```R
# ========================== R code ==========================
imiss.common <- read.table("chr1-22.ASA.maf05.imiss",h=T)
het.common.pruned <- read.table("chr1-22.ASA.maf05.pruned.het",h=T)
imiss.het <- merge(imiss.common, het.common.pruned, by="IID")

up3sd<-mean(imiss.het$F)+3*sd(imiss.het$F)
low3sd<-mean(imiss.het$F)-3*sd(imiss.het$F)
             
plot(imiss.het$F_MISS, imiss.het$F, pch=20, col="darkgrey", xlab="F_MISS",ylab="Inbreeding coefficient (F)")
points(imiss.het$F_MISS[imiss.het$F<low3sd], imiss.het$F[imiss.het$F<low3sd], bg="blue", pch=21)
points(imiss.het$F_MISS[imiss.het$F>up3sd], imiss.het$F[imiss.het$F>up3sd], bg="red", pch=21)
abline(h=up3sd, col="red", lwd=2, lty=2)
abline(h=low3sd, col="blue", lwd=2, lty=2)
legend("bottomleft", c("Below 3SD","Above 3SD"), pt.bg=c("blue","red"), pch=21)
# =============================================================
```
![practical2 het-fmiss](https://user-images.githubusercontent.com/8644480/170840625-0763918e-5686-41ef-b6c5-b3de77468ee6.png)
  
Samples with high missingness usually have abnormal number of heterozygous genotypes called (too many or too few). 
- Remove these samples with high missingness (call rate < 0.98) and redo the heterogeneity check. Use 3 standard deviation (SD) as cut-off and record the samples with abnormally high or low heterogeneity to be removed.

#### **QC: Record the sample IDs with abnormal heterogeneity (outside 3SD) for removal**
```R
# ========================== R code ==========================
imiss.het.mind02 <- imiss.het[imiss.het$F_MISS<0.02,]
up3sd  <- mean(imiss.het.mind02$F)+3*sd(imiss.het.mind02$F)
low3sd <- mean(imiss.het.mind02$F)-3*sd(imiss.het.mind02$F)

excl.up3SD  <- imiss.het.mind02[imiss.het.mind02$F>up3sd,]
excl.low3SD <- imiss.het.mind02[imiss.het.mind02$F<low3sd,]
write.table(rbind(excl.up3SD,excl.low3SD)[,1:2], "to-remove.het3D.indiv", quote=F, row.names=F, col.names=F)
write.table(rbind(excl.up3SD,excl.low3SD)[,2], "to-remove.het3D.iid", quote=F, row.names=F, col.names=F)
# =============================================================
```
  
## Step 4: Duplicated or related individuals
For case-control association analysis, we need to make sure that the samples are biologically unrelated. Usually, samples with second (PI_HAT>0.25) or third degree (PI_HAT>0.125) relatedness are removed.

- Obtain pair-wise IBD for relatedness checking
```bash
~/plink --bfile chr1-22.ASA.maf05.pruned --genome --out chr1-22.ASA.maf05.pruned.IBDcheck
```
- Check the IBD0, IBD1, IBD2 and PI_hat for samples with high missingness
```bash
egrep -wf to-remove.mind02.iid chr1-22.ASA.maf05.pruned.IBDcheck.genome | sort --key 10 -gr | less
```
- Check the IBD0, IBD1, IBD2 and PI_hat for samples with abnormal heterogeneity
```bash
egrep -wf to-remove.het3D.iid chr1-22.ASA.maf05.pruned.IBDcheck.genome | sort --key 10 -gr | less
```
:closed_book: **Q:** Which sample(s) looks like having contamination?
<details>
  <summary> Answer </summary>
  
  - CHSHet002 looks like contaminated sample with high heterogeneity.
</details>

- Check the IBD0, IBD1, IBD2 and PI_hat for other samples
```bash
cat to-remove.*indiv | sort | uniq > to-remove.QC_steps1to3.indiv
cut -d' ' -f 2 to-remove.QC_steps1to3.indiv > to-remove.QC_steps1to3.iid
egrep -wvf to-remove.QC_steps1to3.iid chr1-22.ASA.maf05.pruned.IBDcheck.genome | sort --key 10 -gr | less
```
:closed_book: **Q:** How are the samples in CHSQUAD family related?
<details>
  <summary> Answer </summary>
  
  - CHSQUAD is a quad with F1 being the father and M1 being the mother. C1 and C2 are siblings.
</details>
  
- For those pairs of sample estimated to be closer than second degree kinship (PI_HAT>0.25), remove at least one sample (usually the one with lower call rate) per pair.
<pre>cat > to-remove.related.indiv
CHSQUAD C1
CHSQUAD C2
# the type "Ctrl-D" to quit
</pre>

## Step 5: Population outliers
To validate the self-reported ethnicity and to ensure no population outlier, we merge the genotype data of unrelated samples with the 1000 Genomes Project reference panel and then perform principal component analysis (PCA) using PLINK.
  
- Remove samples with high missingness and related samples and merge with the 1000 Genomes Project data
```bash
cat to-remove.mind02.indiv to-remove.related.indiv > to-remove.mind02_related.indiv
~/plink --bfile chrAll.ASA --remove to-remove.mind02_related.indiv --update-name ASA.1000G.to-update-name.snp --make-bed --out chrAll.ASA.id-1000G.rm-mind02_related
~/plink --bfile chrAll.ASA.id-1000G.rm-mind02_related --bmerge chrAll.ASA.1000GP-All --make-bed --out merged.chrAll.ASA.1000G.rm-mind02_related
```
- Extract only common variants not violating Hardy-Weinberg equilibrium and perform pruning to remove SNPs in LD
```bash
~/plink --bfile merged.chrAll.ASA.1000G.rm-mind02_related \
  --maf 0.05 --hwe 1e-5 --geno 0.05 \
  --indep-pairwise 200 50 0.1 \
  --out merged.chrAll.ASA.1000G.rm-mind02_related
```
- Perform PCA on pruned dataset (PS: You can speed up this step by using more than **1 thread** if your VM is setup with multiple processors)
<pre><code>~/plink --bfile merged.chrAll.ASA.1000G.rm-mind02_related \
  --extract merged.chrAll.ASA.1000G.rm-mind02_related.prune.in \
  --pca 3 \
  --out merged.chrAll.ASA.1000G.rm-mind02_related.pruned \
  <b>--threads 1</b>
</code></pre>
  
- Generate PCA plot
```bash
Rscript practical2.PCAplot.R merged.chrAll.ASA.1000G.rm-mind02_related.pruned.eigenvec
```
<img src="https://user-images.githubusercontent.com/8644480/170877252-273d5367-dbf8-4bd5-8b47-e3cf4b2b17c4.png" width=500>

:closed_book: **Q:** Can you find the population outliers?<br>
:closed_book: **Q:** Can you relate the results of the PCA and heterogeneity analyses for these samples?
<details>
  <summary> Answer </summary>
  
  ```bash
  cut -d' ' -f 2 chrAll.ASA.fam > chrAll.ASA.iid
  awk '$3<0.015' merged.chrAll.ASA.1000G.rm-mind02_related.pruned.eigenvec | egrep -f chrAll.ASA.iid
  cat to-remove.het3D.iid
  ```
</details>  
  
#### **QC: Combine all sample outliers' files and remove these outliers to obtain a new PLINK file with samples passing QC**
```bash
cat to-remove.QC_steps1to3.indiv to-remove.related.indiv > to-remove.QC_steps1to4.indiv
~/plink --bfile chrAll.ASA --remove to-remove.QC_steps1to4.indiv --make-bed --out chrAll.ASA.afterSampleQC
```

- Retry the PCA with removal of all samples not passing QC, including those failed heterogeneity tests in Steps 2 and 3
<details>
  <summary> Try your own PLINK / R codes </summary>

```bash
~/plink --bfile chrAll.ASA.afterSampleQC --update-name ASA.1000G.to-update-name.snp --make-bed --out chrAll.ASA.id-1000G.afterSampleQC
~/plink --bfile chrAll.ASA.id-1000G.afterSampleQC --bmerge chrAll.ASA.1000GP-All --make-bed --out merged.chrAll.ASA.1000G.afterSampleQC
```
<pre><code>
~/plink --bfile merged.chrAll.ASA.1000G.afterSampleQC \
  --maf 0.05 --hwe 1e-5 --geno 0.05 \
  --indep-pairwise 200 50 0.1 \
  --out merged.chrAll.ASA.1000G.afterSampleQC
  
~/plink --bfile merged.chrAll.ASA.1000G.afterSampleQC \
  --extract merged.chrAll.ASA.1000G.afterSampleQC.prune.in \
  --pca 3 \
  --out merged.chrAll.ASA.1000G.afterSampleQC.pruned \
  <b>--threads 1</b>
</code></pre>
  
```bash
Rscript practical2.PCAplot.R merged.chrAll.ASA.1000G.afterSampleQC.pruned.eigenvec 
```
  
<img src="https://user-images.githubusercontent.com/8644480/170875864-98e95f97-b673-4fd9-9126-d4cf545923bd.png" width=500>

</details>

## Variant-QC
It consists of (at least) three steps:

+ Excluding variants with an excessive missing genotype
+ Excluding variants violating Hardy-Weinberg equilibrium (HWE)
+ Excluding variants with a low MAF

The threshold used for filtering depends on sample and genotyping data quality, which vary from study to study. Variant QC should be done carefully as variants removed can be the disease causal variants in which the signals of association may not be completely recovered by imputation.

Here we are using the following thresholds:
+ `--geno 0.02` Call rate <= 98%
+ `--hwe 1e-4`  HWE p < 1x10-4
+ `--maf 0.01`  MAF >= 1%

**QC: We can combine all these variant QCs into one single PLINK command:**
```bash
~/plink --bfile chrAll.ASA.afterSampleQC \
  --geno 0.02 \
  --hwe 1e-4 \
  --maf 0.01 \
  --make-bed \
  --out chrAll.ASA.afterSampleQC.afterVariantQC
```
  
At the end, you will end up with a dataset of 992 samples and 480780 variants passing sample-based and variant-based QCs.

## Association analysis
## Step 1: Examining the phenotype file
Besides the 6th column in PLINK .ped and .fam file, you can use another phenotype file via the command of `--pheno [filename]` in PLINK.

First, let's use R to examine the content of the `CAD_LDL.pheno` phenotype file<br>
:closed_book: **Q1.** How many cornoary artery disease (CAD) patients and controls are there in this dataset? 
<details>
<summary> Answer </summary>

+ 496 cases and 496 controls
```R
# ========================== R code ==========================
pheno <- read.table("CAD_LDL.pheno",h=T)
summary(pheno)
table(pheno$CAD)
# ============================================================
```
</details>
 
:closed_book: **Q2.** Is the low density lipoprotein (LDL) level normally distributed?
<details>
<summary> Answer </summary>

+ Yes. The LDL level is normally distributed and no further transformation is needed
```R
# ========================== R code ==========================
# Plot the histogram to assess if the LDL level is normally distributed
hist(pheno$LDL, freq=T, col="darkred", border ="black", main="Distribution of LDL", xlab="LDL level", ylab="Number of samples")
```
```R
# Plot the QQ plot to assess if the LDL level is normally distributed
qqnorm(pheno$LDL)
qqline(pheno$LDL)
# ============================================================
```
</details>

:closed_book: **Q3.** Is there any relationship between age, LDL level and CAD?
<details>
<summary> Answer </summary>

 + CAD patients are generally older than controls (_P_=2.6x10<sup>-43</sup>) and LDL level increases with age (_P_=0.005).
```R
# ========================== R code ==========================
# Test if age is associated with LDL and CAD
summary(glm(LDL ~ AGE, data=pheno))
summary(glm(CAD==2 ~ AGE, family="binomial", data=pheno))
by(pheno$AGE, pheno$CAD, summary)
# ============================================================
```
</details>

## Step 2: Performing association analysis
Next, we will use PLINK to perform association test
 
### 1. Perform linear regression WITHOUT adjusting for age to test for association with LDL 
```bash  
~/plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
  --chr 1-22,X,XY \
  --pheno CAD_LDL.pheno \
  --pheno-name LDL \
  --linear --adjust \
  --out LDL
```
By default, sex is automatically added as a covariate for all X chromosome SNPs but not anywhere else. The `--adjust` option generates an `.adjusted` file with basic multiple testing corrections (e.g. FDR) for the raw p-values. It also estimates the genomic control factor (lambda) from the data.
 
:closed_book: **Q1.** Which SNP is the top SNP? Can you give the effect size and the confidence interval (CI)?
<details>
  <summary> Answers </summary>
  
+ **Answer 1:** The top SNP is **rs2075650** on chr19 with association _p_-value of 3.8x10<sup>-32</sup>
```bash
head -n 2 LDL.assoc.linear.adjusted 
```
 
+ **Answer 2**: SNP effect beta (95% CI) = -0.9111 (-1.057, -0.7652).<br>
You can add `--ci 0.95` while running the regression test in PLINK or you can compute the 95% CI in R
<pre><code>~/plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
  <b>--snp rs2075650</b> \
  --pheno CAD_LDL.pheno --pheno-name LDL \
  --linear <b>--ci 0.95</b> \
  --out LDL.topsnp 
</code></pre>
 
```R
# ========================== R code ==========================
L95 <- -0.9111-1.96*0.07448
L95
U95 <- -0.9111+1.96*0.07448
U95
# ============================================================
```
</details>  

:closed_book: **Q2.** What is the lamda? Does it indicate inflation of the association statistics?
<details>
  <summary> Answer </summary>

+ lamda=1.00236. No inflation in summary statistics.
</details>
 
- Plot the manhattan plot to visualize the association results  
```bash 
# Rscript practical3.manhattanPlot.R <assoc_file> <outpu_prefix>
Rscript practical3.manhattanPlot.R LDL.assoc.linear manhattanPlot.LDL
```
 
### 2.  Perform linear regression while adjusted for age
To adjust for confounding factor(s), you can use `--covar <filename>` to specify the covariate file, i.e. `CAD_LDL.pheno` in this example. The covariate(s) used for adjustment is specified through `--covar-name`.
 
<pre><code>~/plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
  --chr 1-22,X,XY \
  --pheno CAD_LDL.pheno \
  --pheno-name LDL \
  <b>--covar CAD_LDL.pheno \</b>
  <b>--covar-name AGE \</b>
  <b>--hide-covar \</b>
  --linear --adjust \
  --out LDL.adj-AGE
</code></pre>

- Plot the association results without adjustment of age and compare it with the previous one without adjustment
```bash 
Rscript practical3.manhattanPlot.R LDL.adj-AGE.assoc.linear manhattanPlot.LDL.adj-AGE
```

### 3. Perform logistic regression with and without adjusting for age to test for association with CAD
- Without adjusting for age
```bash 
~/plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
 --chr 1-22,X,XY \
 --pheno CAD_LDL.pheno --pheno-name CAD \
 --logistic \
 --out CAD
```
- Adjusting for age
```bash 
~/plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
 --chr 1-22,X,XY \
 --pheno CAD_LDL.pheno --pheno-name CAD \
 --covar CAD_LDL.pheno --covar-name AGE --hide-covar \
 --logistic \
 --out CAD.adj-AGE
```
- Plot the manhattan plots to compare the difference
```bash 
Rscript practical3.manhattanPlot.R CAD.assoc.logistic manhattanPlot.CAD
Rscript practical3.manhattanPlot.R CAD.adj-AGE.assoc.logistic manhattanPlot.CAD.adj-AGE
```

### Step 3: Annotating the association findings
#### LocusZoom (https://statgen.github.io/localzoom/)<br>
LocalZoom is a tool for generating regional association plot via web browser or offline. For the web-based js version without uploading the summary statistics, a tabix-indexed tab-delimited text file is needed. If you are comfortable uploading your data to a server, you can consider using the my.locuszoom.org upload service with more custom function.
 
Please noted that both the reference (Ref) and alternative (Alt) alleles are needed for visualizing the LD. For simplicity, we use A1 as the alternative allele and A2 as the reference allele here. You will usually flip to the reference strand and get the right Ref/Alt alleles after imputation. 
```bash
# (i) convert "fixed width" PLINK association result format to "tab-delimited"
# (ii) add the A2 allele (specifying the right Ref and Alt alleles are needed for LD plot)
awk 'BEGIN { OFS="\t"; a2["SNP"]="A2" } NR==FNR { a2[$2]=$6 } NR!=FNR { $4=$4"\t"a2[$2]; if (FNR!=1) $2=$1":"$3; print }' \
 chrAll.ASA.afterSampleQC.afterVariantQC.bim \
 LDL.adj-AGE.assoc.linear | \
 bgzip > LDL.adj-AGE.assoc.linear.forLocuszoom.gz
 
# (iii) index the association file
tabix -s 1 -b 3 -e 3 --skip-lines 1 -f LDL.adj-AGE.assoc.linear.forLocuszoom.gz
```
+ First, make sure that `GWAS` and the right genome assembly are chosen, i.e. `GRCh37` for this example dataset
<img src="https://user-images.githubusercontent.com/8644480/171670740-86b20b26-cc20-4512-81cc-d43940be970c.png" width=800>
 
+ Click `Add tab-indexed datafile` > `Browse` to add the LDL.adj-AGE.assoc.linear.forLocuszoom **.gz** and **.gz.tbi** files
<img src="https://user-images.githubusercontent.com/8644480/171672531-00848874-2530-446a-aad4-102475c2ea62.png" width=800>
+ Choose `GWAS scattered plot` as `Type`
<img src="https://user-images.githubusercontent.com/8644480/171672261-3cf5120d-a24b-455a-96e6-959e262b2ef6.png" width=250>
 
+ Under `Variant from columns` > Choose the major allele `A2` as `Ref` allele and the minor allele `A1` as the `Alt` allele
<img src="https://user-images.githubusercontent.com/8644480/171672887-bc1c4203-0880-44a5-9f66-64b376a155a9.png" width=500>

+ Next > Accept option > Type `rs2075650` in the box > Go to the region
+ Choose `LD population: EAS`
+ You can also change the range of region to `19:45300000-45500000` for a zoom in view
<img src="https://user-images.githubusercontent.com/8644480/171678660-f12f0ff9-c4d1-4472-bf9d-089ad58d0164.png" width=800>

#### FUMA (https://fuma.ctglab.nl/)
You can also use the SNP2GENE function in FUMA (Functional Mapping and Annotation of Genome-Wide Association Studies) for annotating, prioritizing, and interpreting the GWAS results. The SNP2GENE function takes GWAS summary statistics as input, and provides functional annotation for all SNPs in genomic areas identified by lead SNPs.
+ SNP2GENE result for [CAD](https://fuma.ctglab.nl/browse/461)
+ SNP2GENE result for [LDL](https://fuma.ctglab.nl/browse/462)
