# hprcPopGenAnalysis

#### Required tools and programming language
- vt (https://genome.sph.umich.edu/wiki/Vt)
- bcftools (http://samtools.github.io/bcftools/bcftools.html)
- plink2 (https://www.cog-genomics.org/plink/2.0/)
- R

#### 0. Download vcfs
```

wget 

```

#### 1. Remove prefix from chromosome names
```
bcftools annotate --rename-chrs hg38.chrNames.txt -O z -o hprc-pggb.grch38.vcf.gz hprc-v1.0-pggb.grch38.1-22+X.vcf.gz

bcftools annotate --rename-chrs chm13.chrNames.txt -O z -o hprc-pggb.chm13.vcf.gz hprc-v1.0-pggb.chm13.1-22+X.vcf.gz

bcftools annotate --rename-chrs hg38.chrNames.txt -O z -o hprc-mc.grch38.vcf.gz hprc-v1.0-mc-grch38.vcf.gz

bcftools annotate --rename-chrs chm13.chrNames.txt -O z -o hprc-mc.chm13.vcf.gz hprc-v1.0-mc-chm13.vcf.gz
```

#### 2. Extract single chromosome from vcf
```
for ass in chm13 hg38; do for c in $(seq 1 22); do 
bcftools view -O z -o hprc-pggb.$ass.chr$c.vcf.gz -r chr$c hprc-pggb.$ass.vcf.gz

bcftools view -O z -o hprc-mc.$ass.chr$c.vcf.gz -r chr$c hprc-mc.$ass.vcf.gz ;
done;done

```

#### 3. Normalize
```
for met in pggb mc; do for c in $(seq 1 22); do
vt normalize -n -r chm13.v1.0.fasta hprc-$met.grch38.chr$c.vcf.gz > hprc-$met.grch38.chr$c.norm.vcf | bgzip hprc-$met.grch38.chr$c.norm.vcf | tabix -p vcf hprc-$met.grch38.chr$c.norm.vcf.gz ;
done; done

for met in pggb mc; do for c in $(seq 1 22); do
vt normalize -n -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hprc-$met.chm13.chr$c.vcf.gz > hprc-$met.chm13.chr$c.norm.vcf | bgzip hprc-$met.chm13.chr$c.norm.vcf | tabix -p vcf hprc-$met.chm13.chr$c.norm.vcf.gz ;
done; done
```

#### 4. Decompose
```
for ass in chm13 hg38; do for c in $(seq 1 22); do
vt decompose hprc-pggb.$met.chr$c.norm.vcf.gz > hprc-pggb.$met.chr$c.deco.vcf | bgzip hprc-pggb.$met.chr$c.deco.vcf | tabix -p vcf hprc-pggb.$met.chr$c.deco.vcf.gz ;
done; done

for ass in chm13 hg38; do for c in $(seq 1 22); do
vt decompose hprc-mc.$met.chr$c.norm.vcf.gz > hprc-mc.$met.chr$c.deco.vcf | bgzip hprc-mc.$met.chr$c.deco.vcf | tabix -p vcf hprc-mc.$met.chr$c.deco.vcf.gz ;
done; done
```

#### 5. Extract SNPs only
```
for ass in chm13 hg38; do for c in $(seq 1 22); do
bcftools view -i 'STRLEN(REF)<=2 & STRLEN(ALT)<=2' -O z -o hprc-pggb.$met.chr$c.SNPs.vcf.gz hprc-pggb.$met.chr$c.deco.vcf.gz | tabix -p vcf hprc-pggb.$met.chr$c.SNPs.vcf.gz ;
done;done

for ass in chm13 hg38; do for c in $(seq 1 22); do
bcftools view -i 'STRLEN(REF)<=2 & STRLEN(ALT)<=2' -O z -o hprc-mc.$met.chr$c.SNPs.vcf.gz hprc-mc.$met.chr$c.deco.vcf.gz |tabix -p vcf hprc-mc.$met.chr$c.SNPs.vcf.gz;
done;done
```
#### 6. Remove reference from samples
```
for met in pggb mc; do for c in $(seq 1 22); do
vcftools --gzvcf hprc-$met.chm13.chr$c.SNPs.vcf.gz --remove-indv grch38 --out hprc-$met.chm13.chr$c.SNPs.noRef --recode --keep-INFO-all | bgzip hprc-$met.chm13.chr$c.SNPs.noRef.recode.vcf | tabix -p vcf hprc-$met.chm13.chr$c.SNPs.noRef.recode.vcf.gz; 
done;done

for met in pggb mc; do for c in $(seq 1 22); do
vcftools --gzvcf hprc-$met.grch38.chr$c.SNPs.vcf.gz --remove-indv chm13 --out hprc-$met.grch38.chr$c.SNPs.noRef --recode --keep-INFO-all | bgzip hprc-$met.grch38.chr$c.SNPs.noRef.recode.vcf | tabix -p vcf hprc-$met.grch38.chr$c.SNPs.noRef.recode.vcf.gz ;
done;done
```
#### 7. PCA for entire chromosome
```
for ass in chm13 hg38; do for c in $(seq 1 22); do
plink2 --vcf hprc-pggb.$ass.chr$c.SNPs.noRef.recode.vcf.gz --double-id --set-all-var-ids @:#$r:$a --rm-dup exclude-mismatch --vcf-half-call m --maf 0.01 --freq --out hprc-pggb.$ass.chr$c.SNPs | plink2 --vcf hprc-pggb.$ass.chr$c.SNPs.noRef.recode.vcf.gz --double-id --set-all-var-ids @:#$r:$a --rm-dup exclude-mismatch --vcf-half-call m --make-bed --read-freq hprc-pggb.$ass.chr$c.SNPs.afreq --pca --out hprc-pggb.$ass.chr$c.SNPs ;
done; done

for ass in chm13 hg38; do for c in $(seq 1 22); do
plink2 --vcf hprc-mc.$ass.chr$c.SNPs.noRef.recode.vcf.gz --double-id --set-all-var-ids @:#$r:$a --rm-dup exclude-mismatch --vcf-half-call m --maf 0.01 --freq --out hprc-mc.$ass.chr$c.SNPs | plink2 --vcf hprc-mc.$ass.chr$c.SNPs.noRef.recode.vcf.gz --double-id --set-all-var-ids @:#$r:$a --rm-dup exclude-mismatch --vcf-half-call m --make-bed --read-freq hprc-mc.$ass.chr$c.SNPs.afreq --pca --out hprc-mc.$ass.chr$c.SNPs ;
done; done
```
#### 8. Extract variants on p arm
```
for ass in chm13 hg38; do for c in $(seq 1 22); do
bcftools view -R chr$c.$ass.pArm_coord.txt -O z -o hprc-pggb.$ass.chr$c.SNPs.pArm.vcf.gz hprc-pggb.$ass.chr$c.SNPs.noRef.recode.vcf.gz | tabix -p vcf hprc-pggb.$ass.chr$c.SNPs.pArm.vcf.gz;
done; done

for ass in chm13 hg38; do for c in $(seq 1 22); do
bcftools view -R chr$c.$ass.pArm_coord.txt -O z -o hprc-mc.$ass.chr$c.SNPs.pArm.vcf.gz hprc-mc.$ass.chr$c.SNPs.noRef.recode.vcf.gz | tabix -p vcf hprc-mc.$ass.chr$c.SNPs.pArm.vcf.gz;
done; done
```
#### 9. Extract variants on q arm
```
for ass in chm13 hg38; do for c in $(seq 1 22); do
bcftools view -R chr$c.$ass.qArm_coord.txt -O z -o hprc-pggb.$ass.chr$c.SNPs.qArm.vcf.gz hprc-pggb.$ass.chr$c.SNPs.noRef.recode.vcf.gz | tabix -p vcf hprc-pggb.$ass.chr$c.SNPs.qArm.vcf.gz;
done; done

for ass in chm13 hg38; do for c in $(seq 1 22); do
bcftools view -R chr$c.$ass.qArm_coord.txt -O z -o hprc-mc.$ass.chr$c.SNPs.qArm.vcf.gz hprc-mc.$ass.chr$c.SNPs.noRef.recode.vcf.gz | tabix -p vcf hprc-mc.$ass.chr$c.SNPs.qArm.vcf.gz;
done; done
```
#### 10. PCA using variants on p and q arms

#### 11. Cluster analysis
```
Rscript scr/clusterAnalysis.R pathToImputFiles length assembly method pathToOutputFile
```

#### 12. PCA plot
