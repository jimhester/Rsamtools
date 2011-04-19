TabixFile <-
    function(file, index=paste(file, "tbi", sep="."), ...)
{
    .RsamtoolsFile(.TabixFile, file, index)
}

open.TabixFile <-
    function(con, ...)
{
    .io_check_exists(index(con))
    ## FIXME: path? index?
    con$.extptr <- .Call(.tabixfile_open, path(con), index(con))
    invisible(con)
}

close.TabixFile <-
    function(con, ...)
{
    if (!isOpen(con))
        stop("isOpen(<TabixFile>) is not 'TRUE'")
    con$.extptr <- .Call(.tabixfile_close, .extptr(con))
    invisible(con)
}

setMethod(isOpen, "TabixFile",
    function(con, rw="")
{
    if (!missing(rw) && rw != "read")
        stop("'rw' must be 'read'")
    .Call(.tabixfile_isopen, .extptr(con))
})

bgzipTabix <-
    function(fromFname, toFname = paste(fromFname, "gz", sep="."),
             overwrite=FALSE)
{
    .Call(.bgzip_tabix, fromFname, toFname, overwrite)
}

indexTabix <- 
    function(file,
             format=c("gff", "bed", "sam", "vcf", "vcf4", "psltbl"),
             seq=integer(), begin=integer(), end=integer(),
             skip=0L, comment="#", zeroBased=FALSE, ...)
{
    tryCatch({
        format <- 
            if (!missing(format)) match.arg(format)
            else character()
        idx <- .Call(.index_tabix, file, format,
                     seq, begin, end, skip, comment, zeroBased)
        sprintf("%s.tbi", file)
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", file)
    })
}

.seqnamesTabix <-
    function(file, ...)
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    .Call(.seqnames_tabix, .extptr(file))
}

setMethod(seqnamesTabix, "TabixFile", .seqnamesTabix)

setMethod(seqnamesTabix, "character", function(file, ...) {
    .seqnamesTabix(TabixFile(file))
})

.tabix_scan <-
    function(file, ..., space, start, end)
{
    tryCatch({
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }

        yieldSize <- 1000000L           # a guess, grows as necessary
        result <- .Call(.scan_tabix, .extptr(file),
                        list(space, start, end), yieldSize)
        names(result) <- sprintf("%s:%d-%d", space, start, end)
        result
    }, error=function(err) {
        stop("scanTabix: ", conditionMessage(err), "\n  path: ",
             path(file), call.=FALSE)
    })
}

setMethod(scanTabix, c("TabixFile", "RangesList"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., space=space(param),
                start=.uunlist(start(param)),
                end=.uunlist(end(param)))
})

setMethod(scanTabix, c("TabixFile", "RangedData"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., param=ranges(param))
})

setMethod(scanTabix, c("TabixFile", "GRanges"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., space=as.character(seqnames(param)),
                start=start(param), end=end(param))
})

setMethod(scanTabix, c("character", "ANY"),
    function(file, ..., param)
{
    file <- TabixFile(file)
    scanTabix(file, ..., param=param)
})

.tabix_yield <-
    function(file, ..., yieldSize)
{
    tryCatch({
        if (!isOpen(file))
            open(file)
        .Call(.yield_tabix, .extptr(file), as.integer(yieldSize))
    }, error=function(err) {
        stop("yield: ", conditionMessage(err), "\n  path: ",
             path(file), call.=FALSE)
    })
}

setMethod(yieldTabix, "TabixFile",
    function(file, ..., yieldSize=1000000L)
{
    .tabix_yield(file, ..., yieldSize=yieldSize)
})