destfile <- tempfile()
download.file("https://github.com/privefl/bigsnpr/raw/master/tmp-tests/reprod.rds",
              destfile = destfile)
DATA <- readRDS(destfile)
str(DATA)

X <- DATA$data
table(X, exclude = NULL)
y <- DATA$label
table(y, exclude = NULL)

require(xgboost)
bst <- xgboost(data = X, label = y, nrounds = 20, save_period = NULL)
dt <- xgb.model.dt.tree(model = bst)
table(dt$Split)

set.seed(1)
X.NA <- X
len <- length(X.NA)
X.NA[sample(len, 0.1*len)] <- NA
table(X.NA, exclude = NULL)

bstNA <- xgboost(data = X.NA, label = y, nrounds = 20, save_period = NULL)
dtNA <- xgb.model.dt.tree(model = bstNA)
table(dtNA$Split)
