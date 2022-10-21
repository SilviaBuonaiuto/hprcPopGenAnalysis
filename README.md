# hprcPopGenAnalysis

#### Required tools and programming language
- vt (https://genome.sph.umich.edu/wiki/Vt)
- bcftools (http://samtools.github.io/bcftools/bcftools.html)
- plink2 (https://www.cog-genomics.org/plink/2.0/)
- R

#### 0. Download vcfs
```
for c in $(seq 1 22); do 
wget http://hypervolu.me/~guarracino/HPRC/wgg.88_confident_variants_APR_08_2022/chr1.pan.fa.a2fb268.4030258.6a1ecc2.smooth.reliable.vcf.gz* ; \
 done

```

#### 1. Remove prefix from chromosome names
```
for c in $(seq 1 22); do 
bcftools annotate --rename-chrs chr$c.hg38.chrNames.txt \
-O z \
-o chr$c.confident.vcf.gz \
chr$c.pan.fa.a2fb268.4030258.6a1ecc2.smooth.reliable.vcf.gz; done
```

#### 2. Normalize
```
for c in $(seq 1 22); do
vt normalize \
-n \ 
-r chm13.v1.0.fasta chr$c.confident.vcf.gz > chr$c.confident.norm.vcf | \
bgzip chr$c.confident.norm.vcf | \
tabix -p vcf chr$c.confident.norm.vcf.gz ; done
```

#### 3. Decompose
```
for c in $(seq 1 22); do 
vt decompose chr$c.confident.norm.vcf.gz > chr$c.confident.deco.vcf | \
bgzip chr$c.confident.deco.vcf | \
tabix -p vcf chr$c.confident.deco.vcf.gz ; done
```

#### 4. Extract SNPs only
```
for c in $(seq 1 22); do
bcftools view -i 'STRLEN(REF)<=2 & STRLEN(ALT)<=2' \
-O z \
-o chr$c.confident.SNPs.vcf.gz chr$c.confident.deco.vcf.gz | \
tabix -p vcf chr$c.confident.SNPs.vcf.gz ; done
```
#### 5. Remove reference from samples
```
for c in $(seq 1 22); do
vcftools --gzvcf chr$c.confident.SNPs.vcf.gz \
--remove-indv grch38 \
--out chr$c.confident.SNPs.noRef \
--recode \
--keep-INFO-all | \
bgzip chr$c.confident.SNPs.noRef.recode.vcf | \
tabix -p vcf chr$c.confident.SNPs.noRef.recode.vcf.gz; done
```
#### 6. PCA for entire chromosome
```
for ass in chm13 hg38; do for c in $(seq 1 22); do 
plink2 --vcf chr$c.confident.SNPs.noRef.recode.vcf.gz\
--double-id \
--set-all-var-ids @:#$r:$a \
--rm-dup exclude-mismatch \
--vcf-half-call m \
--maf 0.01 \
--freq \
--out chr$c.confident.SNPs | \
plink2 --vcf chr$c.confident.SNPs.noRef.recode.vcf.gz \
--double-id \
--set-all-var-ids @:#$r:$a \
--rm-dup exclude-mismatch \
--vcf-half-call m \
--make-bed \
--read-freq chr$c.confident.SNPs.afreq \
--pca \
--out chr$c.confident.SNPs ; done
```
#### 7. Extract variants on p arm
```
for c in $(seq 1 22); do 
bcftools view -R chr$c.chm13.pArm_coord.txt \
-O z \
-o chr$c.confident.SNPs.pArm.vcf.gz chr$c.confident.SNPs.noRef.recode.vcf.gz | \
tabix -p vcf chr$c.confident.SNPs.pArm.vcf.gz; done
```
#### 8. Extract variants on q arm
```
for ass in chm13 hg38; do for c in $(seq 1 22); do
bcftools view -R chr$c.chm13.qArm_coord.txt \
-O z \
-o chr$c.confident.SNPs.qArm.vcf.gz \
chr$c.confident.SNPs.noRef.recode.vcf.gz | \
tabix -p vcf chr$c.confident.SNPs.qArm.vcf.gz; done
```
#### 9. PCA using variants on p and q arms

#### 10. Cluster analysis and plot
```
Rscript scr/clusterAnalysis.R \
pathToImputFiles \
length \
assembly \
method \
pathToOutputFile


Rscript scr/plotClustersHprc.R \
pathToImputFiles \
pathToPlot
```

#### 11. PCA plot
```
## put together eigenvec files
python3 scr/appendEigenvec.py \
-i "pathToImputFiles" \
-met method \
-ass assembly \
-len lehgth \
-o pathToOutputFile

Rscript scr/PCA_plot.R "pathToImputFiles" \
met \
pathToOtputPlot
```