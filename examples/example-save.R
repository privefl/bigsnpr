set.seed(1)

# Reading example
test <- snp_attachExtdata()

# I can add whatever I want to an S3 class
test$map$`p-values` <- runif(nrow(test$map))
str(test$map)

# Reading again
test.savedIn <- sub("\\.bk$", ".rds", test$genotypes$backingfile)
test2 <- snp_attach(rdsfile = test.savedIn)
str(test2$map) # new slot wasn't saved

# Save it
test <- snp_save(test)

# Reading again
test3 <- snp_attach(rdsfile = test.savedIn)
str(test3$map) # it is saved now

# The complicated code of this function
snp_save
