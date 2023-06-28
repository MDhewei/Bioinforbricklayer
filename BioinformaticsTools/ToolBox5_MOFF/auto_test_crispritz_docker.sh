#!/bin/sh

#CREATE TEST DIRECTORY
mkdir test_crispritz
cd test_crispritz/

echo "DOWNLOADING ANNOTATIONS FILES"
#ANNOTATIONS FILE PATH AND ANNOTATIONS BED DIRECTORY DOWNLOAD
curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/hg38Annotation.zip --output hg38Annotation.zip
unzip hg38Annotation.zip
echo "ANNOTATIONS FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING GUIDES FILES"
#GUIDE DIRECTORY DOWNLOAD
curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/EMX1.sgRNA.txt --output EMX1.sgRNA.txt
echo "GUIDES FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING PAM FILES"
#PAM DIRECTORY DOWNLOAD
curl https://raw.githubusercontent.com/pinellolab/CRISPRitz/master/test_scripts/20bp-NGG-SpCas9.txt --output 20bp-NGG-SpCas9.txt
echo "PAM FILES DOWNLOADED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING CHR22 (TIME DEPENDING ON YOUR INTERNET CONNECTION)"
#CHR22 DOWNLOAD AND gunzip --quiet -q
mkdir hg38_ref
cd hg38_ref/
curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz --output chr22.fa.gz
gunzip --quiet chr22.fa.gz
cd ..
echo "CHR22 DOWNLOADED AND EXTRACTED"
echo "-------------------------------------------------------------------------"

echo "DOWNLOADING CHR22_VCF (TIME DEPENDING ON YOUR INTERNET CONNECTION)"
#CHR22 VCF DOWNLOAD
mkdir hg38_1000genomeproject_vcf
cd hg38_1000genomeproject_vcf/
curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz --output ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
cd ..
echo "CHR22_VCF DOWNLOADED AND EXTRACTED"
echo "-------------------------------------------------------------------------"

echo "EVERYTHING DOWNLOADED AND READY TO USE"

echo ""
echo "STARTING THE AUTO TEST SCRIPT"
echo "THIS SCRIPT WILL NOW TEST ALL THE CRISPRitz FUNCTION TO CHECK THE INSTALLATION"
echo ""
echo "TESTING ADD-VARIANTS"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py add-variants hg38_1000genomeproject_vcf/ hg38_ref/ &>output.redirect.out
echo -e "ADD-VARIANTS \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING INDEX-GENOME"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py index-genome hg38_ref hg38_ref/ 20bp-NGG-SpCas9.txt -bMax 2 &>output.redirect.out
echo -e "INDEX-GENOME \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING SEARCH WITH ONLY MISMATCHES"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py search hg38_ref/ 20bp-NGG-SpCas9.txt EMX1.sgRNA.txt emx1.hg38 -mm 4 -t -scores hg38_ref/ &>output.redirect.out
echo -e "SEARCH WITH ONLY MISMATCHES \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING SEARCH WITH MISMATCHES AND BULGES"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py search genome_library/NGG_2_hg38_ref/ 20bp-NGG-SpCas9.txt EMX1.sgRNA.txt emx1.hg38.bulges -index -mm 4 -bDNA 1 -bRNA 1 -t &>output.redirect.out
echo -e "SEARCH WITH MISMATCHES AND BULGES \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING ANNOTATE-RESULTS"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py annotate-results emx1.hg38.targets.txt hg38Annotation.bed emx1.hg38 &>output.redirect.out
echo -e "ANNOTATE-RESULTS \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

echo "TESTING GENERATE-REPORT"
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispritz crispritz.py generate-report GAGTCCGAGCAGAAGAAGAANNN -mm 4 -annotation emx1.hg38.Annotation.summary.txt -extprofile emx1.hg38.extended_profile.xls &>output.redirect.out
echo -e "GENERATE-REPORT \e[32mOK\e[0m"
echo "-------------------------------------------------------------------------"

rm output.redirect.out
echo ""
echo -e "EVERY TEST \e[32mPASSED\e[0m!!! ENJOY CRISPRitz"
