.filter_sort_mtx_files <- function(namesvec) {
    files <- .FEATURE_BC_MATRIX_FILES
    names(files) <- files
    res <- lapply(files, function(file) {
        namesvec[startsWith(namesvec, file)]
    })
    unlist(res)
}

.check_filter_mtx <- function(filelist) {
    afiles <- .filter_sort_mtx_files(names(filelist))
    if (!identical(names(afiles), .FEATURE_BC_MATRIX_FILES))
        stop(
            "'TENxFileList' does not contain the expected files:\n  ",
            .FEATURE_BC_MATRIX_FILES_PRINT
        )
    filelist[afiles]
}

.file_for_format <- function(fdir, format, ...) {
    if (identical(format, "h5")) {
        h5f <- file.path(fdir, paste0(.FEATURE_MATRIX_FILE_STUB, ".", format))
        if (!file.exists(h5f))
            stop("The '", basename(h5f), "' file was not found.")
        path <- TENxH5(h5f, ...)
    } else if (identical(format, "mtx")) {
        mtxf <- file.path(
            fdir, paste0(.FEATURE_MATRIX_FILE_STUB, ".tar.gz")
        )
        if (!file.exists(mtxf))
            stop("The '", basename(mtxf), "' file was not found.")
        path <- TENxFileList(mtxf)
    }
    path
}

.filter_h5_files <- function(path, format) {
    fname <- file.path(path, paste0(.FEATURE_MATRIX_FILE_STUB, ".", format))
    ish5file <- endsWith(names(path), fname)
    if (!any(ish5file))
        stop("The '", fname, "' file was not found.")
    path(path[ish5file])
}

#' @importFrom TENxIO TENxFileList TENxH5
.find_convert_resources <- function(path, format, ...) {
    if (!is(path, "TENxFileList"))
        path <- .file_for_format(path, format, ...)
    if (identical(format, "h5"))
        path <- .filter_h5_files(path, format)
    else if (identical(format, "mtx"))
        path <- .check_filter_mtx(path)

    path
}

#' @importFrom BiocIO FileForFormat
.filter_xenium_file <- function(path) {
    xf <- list.files(path, pattern = "\\.xenium$", full.names = TRUE)
    FileForFormat(xf)
}

#' @importClassesFrom VisiumIO TENxSpatialParquet
.filter_parquet_file <- function(path) {
    pf <- list.files(path, pattern = "cells\\.parquet$", full.names = TRUE)
    FileForFormat(pf, prefix = "TENxSpatial", suffix = NULL)
}

.boundaries_for_format <- function(fdir, fileext) {
    fname <- paste0("cell_boundaries", ".", fileext)
    cellsf <- file.path(fdir, fname)
    if (!file.exists(cellsf))
        stop("The '", basename(cellsf), "' file was not found.")
    ## override format for csv.gz files
    format <-
        if (identical(fileext, "csv.gz")) "csv" else tools::file_ext(cellsf)
    FileForFormat(
        cellsf, format = format, prefix = "TENxSpatial", suffix = NULL
    )
}

#' @importFrom BiocBaseUtils checkInstalled
.cache_url_file <- function(url, redownload = FALSE) {
    checkInstalled("BiocFileCache")
    bfc <- BiocFileCache::BiocFileCache()
    bquery <- BiocFileCache::bfcquery(bfc, url, "rname", exact = TRUE)
    ## only re-download manually b/c bfcneedsupdate always returns TRUE
    if (identical(nrow(bquery), 1L) && redownload)
        BiocFileCache::bfcdownload(
            x = bfc, rid = bquery[["rid"]], ask = FALSE
        )

    BiocFileCache::bfcrpath(
        bfc, rnames = url, exact = TRUE, download = TRUE, rtype = "web"
    )
}
