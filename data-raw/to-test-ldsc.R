#### Run bulik/ldsc ####
docker run -v "/${PWD}":/home/ldsc_local --rm -ti manninglab/ldsc

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/PASS_BMI1.sumstats.gz

ldsc.py \
--h2 PASS_BMI1.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ../ldsc_local/bmi
# Total Observed scale h2: 0.1408 (0.0066)
# Intercept: 0.772 (0.0082)

ldsc.py \
--h2 PASS_BMI1.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ../ldsc_local/bmi \
--no-intercept
# Total Observed scale h2: 0.039 (0.0057)

ldsc.py \
--h2 PASS_BMI1.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ../ldsc_local/bmi \
--two-step 10
# Total Observed scale h2: 0.1381 (0.0061)
# Intercept: 0.7772 (0.0069)

#### Run in R ####
library(bigreadr)
sumstats <- fread2("https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/PASS_BMI1.sumstats.gz")

M <- sum(fread2(paste0("eur_w_ld_chr/", 1:22, ".l2.M_5_50")))
ld <- fread2(paste0("eur_w_ld_chr/", 1:22, ".l2.ldscore.gz"))

sumstats$L2 <- ld$L2[match(sumstats$SNP, ld$SNP)]
str(sumstats2 <- na.omit(sumstats))

input <- list(ld_score = sumstats2$L2,
              ld_size = M,
              chi2 = sumstats2$CHISQ,
              sample_size = sumstats2$N)
# saveRDS(input, "tmp-data/sumstats-ldsc.rds")

input$ncores <- 4
do.call(bigsnpr::snp_ldsc, args = input)
input$intercept <- 1
do.call(bigsnpr::snp_ldsc, args = input)
