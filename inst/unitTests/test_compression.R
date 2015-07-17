test_bgzip_openclose <- function()
{
    ## trying to determine that file handle has been cleaned up
    checkIdentical(TRUE, dir.create(d <- tempfile()))
    fin <- file.path(d, "in")
    fout <- file.path(d, "out")
    writeLines("123", con=fin)
    bgzip(fin, fout)
    checkIdentical(TRUE, file.remove(fin))
    checkIdentical(TRUE, file.remove(fout))
    checkIdentical(0L, unlink(d, recursive=TRUE))
}

## NH 7/16/2015: I don't know what this is about, so I'm leaving it alone
test_razip_small_files <- function()
{
    src <- system.file("extdata", "ce2dict1.fa", package="Rsamtools")
    file.copy(src, dest <- tempfile())
    checkIdentical(readLines(src), readLines(dest))
}
