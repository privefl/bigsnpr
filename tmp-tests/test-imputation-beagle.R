# Convert bed/bim/fam to vcf
# ./plink -recode vcf -bfile ../../Bureau/POPRES_data/popresSub

# Impute with BEAGLE 4.1
# java -Xmx6g -jar ../beagle.21Jan17.6cc.jar gt=plink.vcf out=test
