set.seed(1)

#reading example
test <- snp_attachExtdata()

# I can add whatever I want to an S3 class
test$map$`p-values` <- runif(nrow(test$map))
str(test$map)

# reading again
test2 <- snp_attach(rdsfile = test$savedIn)
str(test2$map) # new slot wasn't saved

# save it
test2$map$`p-values` <- runif(nrow(test2$map))
test2 <- snp_saveModifs(test2)

# reading again
test3 <- snp_attach(rdsfile = test2$savedIn)
str(test3$map) # it is saved now

# the complicated code of this function
snp_saveModifs
