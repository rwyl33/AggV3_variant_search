# AggV3_variant_search
This script finds GCGR variants inherited as alternate homozygous within AggV3, joins genotype information with information on functional annotation and participant information. 

```
#This intersects coordinates of GCGR gene body with BED file which contains filepaths to shards. 
#This will output the shard and subshard that contains variants in the GCGR gene
bedtools intersect -wo -a GCGR_gene_body.bed -b filesystems/biallelic_shards.bed
#This finds the GCGR gene body to be located in shard 82, subshard 14 and gives the s3 link to its location

#The subshard containing the GCGR gene body is mounted into the session
#This filters for all samples within GCGR gene body which PASS QC and have at least one sample which is alternate homozygous
bcftools view -r chr17:81804150-81814008 -f PASS \
  -S filesystems/aggv3_consented_samples.txt \
  filesystems/dragen.vcf.gz -Oz -o consented_subset.vcf.gz

#This filters consenting samples for variants which are alternate homozygous
bcftools view -i 'GT[*]="1/1"' consented_subset.vcf.gz -Oz -o filtered_genotypes.vcf.gz

#This filters samples for those which have homozygous alternate variants and puts them into a tsv
bcftools query -HH -f '[%SAMPLE\t%CHROM:%POS:%REF-%ALT\t%FILTER\t%GT\n]' \
  filtered_genotypes.vcf.gz | awk -F'\t' 'NR==1 || $4 == "1/1"' > hom_alt.tsv
wc -l hom_alt.tsv 

#Take gene annotations in the region chr17:81804150-81814008 and filter for those which are gnomADg_AF_genomes < 0.02 
bedtools intersect -wo -a GCGR_gene_body.bed -b filesystems/functional_annotation_shards.bed

bcftools +split-vep filesystems/dragen.gel.annotated.vcf.gz -l

#Extract variants in GCGR
bcftools +split-vep filesystems/dragen.gel.annotated.vcf.gz \
    -HH \
    -r chr17:81804150-81814008 \
    -f '%CHROM:%POS:%REF-%ALT\t%EXON\t%INTRON\t%Protein_position\t%Amino_acids\t%gnomADg_AF_genomes\t%CADD_PHRED\t%CADD_RAW\t%ClinVar\t%PhyloP\t%SpliceAI_pred_DS_AG\t%SpliceAI_pred_DS_AL\t%SpliceAI_pred_DS_DG\t%SpliceAI_pred_DS_DL\n' \
    -d \
    > GCGR_annotations.tsv

#filter annotated variants for those which are gnomADg_AF_genomes < 0.02
awk -F'\t' 'NR==1 || ($6!="." && $6<0.02)' GCGR_annotations.tsv > GCGR_annotations_rare.tsv
wc -l GCGR_annotations_rare.tsv

#sort both genotype and annotated variant files
sort -k2,2 hom_alt.tsv > hom_alt_sorted.tsv
sort -k1,1 GCGR_annotations_rare.tsv > annotations_rare_sorted.tsv

#join list of homozygous alternate variants (genotype file) with annotated rare variants (AF < 0.02)
join -t $'\t' -1 2 -2 1 hom_alt_sorted.tsv annotations_rare_sorted.tsv > joined.tsv
wc -l joined.tsv

#Add headers back to joined.tsv
sed -i '1i VARIANT\tSAMPLE\tFILTER\tGT\tEXON\tINTRON\tProtein_position\tAmino_acids\tgnomADg_AF_genomes\tCADD_PHRED\tCADD_RAW\tClinVar\tPhyloP\tSpliceAI_pred_DS_AG\tSpliceAI_pred_DS_AL\tSpliceAI_pred_DS_DG\tSpliceAI_pred_DS_DL' joined.tsv

#Combine genotype-annotation file with participant information
csvformat -t joined.tsv > joined_as_csv.csv

conda install csvkit
csvjoin -I -c "SAMPLE,platekey" joined_as_csv.csv filesystems/sample_list_aggv3.csv \
  | csvcut -c "SAMPLE,VARIANT,FILTER,GT,participant_id,type,study_source,EXON,INTRON,Protein_position,Amino_acids,gnomADg_AF_genomes,CADD_PHRED,CADD_RAW,ClinVar,PhyloP,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL" \
  | csvformat -T > combined.tsv

#Extract unique variants (given in column 2) of combined.tsv into .txt file for airlock
awk -F'\t' 'NR>1 {print $2}' combined.tsv | sort -u > variants_airlock.txt
wc -l variants_airlock.txt

#This prints all unique participant IDs into a txt file (required for airlock)
awk -F'\t' 'NR>1 {print $5}' combined.tsv | sort -u > participant_ids_airlock.txt
wc -l participant_ids_airlock.txt

#This deduplicates variants and puts into a tsv 
awk -F'\t' '!seen[$2]++' combined.tsv > variants.tsv

awk -F'\t' '$12!="." && $12!="" && $12!="-" && $12<0.02' variants.tsv | wc -l
awk -F'\t' '$12!="." && $12!="" && $12<0.01' variants.tsv | wc -l
awk -F'\t' '$12!="." && $12!="" && $12<0.001' variants.tsv | wc -l
awk -F'\t' '$12!="." && $12!="" && $12<0.0001' variants.tsv | wc -l
awk -F'\t' '$12!="." && $12!="" && $12<0.00001' variants.tsv | wc -l
awk -F'\t' '$12!="." && $12!="" && $12<0.000001' variants.tsv | wc -l

#use variants.tsv to filter CADD_PHRED score
head -1 variants.tsv | tr '\t' '\n' | cat -n
#CADD_PHRED is given in column 13
(head -1 variants.tsv; tail -n +2 variants.tsv \
  | sort -t$'\t' -k13,13 -nr variants.tsv | head -20 ) > CADD_PHRED_sorted.tsv
#AggV3 annotated using CADD v1.6 plugin

```
