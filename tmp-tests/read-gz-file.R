write.table(data.frame(A1 = letters,
                       ID1 = round(runif(26), 3),
                       ID2 = round(runif(26), 3))[rep(1:26, 10), ],
            gzfile("tmp-data/test.txt.gz"), row.names = FALSE)

readBin("tmp-data/test.txt.gz", what = raw(), n = 10)

# con <- gzfile("tmp-data/test.txt.gz")
# str(con)
# object.size(con)

bigreadr::nlines("tmp-data/test.txt.gz")
vroom::vroom("tmp-data/test.txt.gz")
