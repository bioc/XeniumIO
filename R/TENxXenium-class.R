#' @include XeniumFile.R
#' @importClassesFrom TENxIO TENxFileList TENxH5
setClassUnion("TENxFileList_OR_TENxH5", members = c("TENxFileList", "TENxH5"))

setClassUnion(
    "TENxSpatialParquet_OR_TENxSpatialCSV",
    members = c("TENxSpatialParquet", "TENxSpatialCSV")
)

#' @docType class
#'
#' @title A class to represent Xenium output data
#'
#' @description This class is a composed class of
#'   [TENxFileList][TENxIO::TENxFileList-class] which can contain a list of
#'   [TENxFile][TENxIO::TENxFile-class] objects for the cell-feature matrix. It
#'   is meant to handle a single Xenium sample from 10X Genomics.
#'
#' @slot resources A [TENxFileList][TENxIO::TENxFileList-class] or
#'   [TENxH5][TENxIO::TENxH5] object containing the cell feature matrix.
#'
#' @slot boundaries Either a
#'   [TENxSpatialParquet][VisiumIO::TENxSpatialParquet-class] or
#'   [TENxSpatialCSV][VisiumIO::TENxSpatialCSV-class] object containing the
#'   spatial boundaries data.
#'
#' @slot coordNames `character()` A vector specifying the names
#'   of the columns in the spatial data containing the spatial coordinates.
#'
#' @slot sampleId `character(1)` A scalar specifying the sample identifier.
#'
#' @slot colData `TENxSpatialParquet` A
#'   [TENxSpatialParquet][VisiumIO::TENxSpatialParquet-class] object containing
#'   the spatial coordinates data.
#'
#' @slot metadata `XeniumFile` A [XeniumFile][XeniumIO::XeniumFile-class] object
#'   containing the metadata information.
#'
#' @return A [SpatialExperiment][SpatialExperiment::SpatialExperiment-class]
#'   object
#'
#' @seealso <https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs>
#'
#' @examples
#' showClass("TENxXenium")
#'
#' @exportClass TENxXenium
.TENxXenium <- setClass(
    Class = "TENxXenium",
    slots = c(
        resources = "TENxFileList_OR_TENxH5",
        boundaries = "TENxSpatialParquet_OR_TENxSpatialCSV",
        coordNames = "character",
        sampleId = "character",
        colData = "TENxSpatialParquet",
        metadata = "XeniumFile"
    )
)

.FEATURE_BC_MATRIX_FILES <- c("barcodes.tsv", "features.tsv", "matrix.mtx")

.FEATURE_BC_MATRIX_FILES_PRINT <- paste0(
    sQuote(.FEATURE_BC_MATRIX_FILES), collapse = ", "
)

.FEATURE_MATRIX_FILE_STUB <- "cell_feature_matrix"

#' @rdname TENxXenium-class
#'
#' @inheritParams VisiumIO::TENxVisium
#'
#' @param xeniumOut `character(1)` The path to the Xenium output directory.
#'
#' @param boundaries_format `character(1)` Either "parquet" or "csv.gz" to
#'   specify the file extension of the boundaries file. Default is "parquet".
#'
#' @details Note that one can provide a `ref` argument to `import` method which
#'   will get passed to the internal `splitAltExps` operation. This allows one
#'   to set a `mainExpName` in the output object.
#'
#' @importFrom methods is new
#' @importFrom BiocBaseUtils isScalarCharacter
#'
#' @examples
#' zipfile <- paste0(
#'     "https://mghp.osn.xsede.org/bir190004-bucket01/BiocXenDemo/",
#'     "Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip"
#' )
#' destfile <- XeniumIO:::.cache_url_file(zipfile)
#' outfold <- file.path(
#'     tempdir(), tools::file_path_sans_ext(basename(zipfile))
#' )
#' if (!dir.exists(outfold))
#'     dir.create(outfold, recursive = TRUE)
#' unzip(
#'     zipfile = destfile, exdir = outfold, overwrite = FALSE
#' )
#' TENxXenium(xeniumOut = outfold) |>
#'     import(ref = "Gene Expression")
#' @export
TENxXenium <- function(
    resources,
    xeniumOut,
    sample_id = "sample01",
    format = c("mtx", "h5"),
    boundaries_format = c("parquet", "csv.gz"),
    spatialCoordsNames = c("x_centroid", "y_centroid"),
    ...
) {
    format <- match.arg(format)
    boundaries_format <- match.arg(boundaries_format)

    if (!missing(xeniumOut)) {
        if (isScalarCharacter(xeniumOut) && !dir.exists(xeniumOut))
            stop(
                "The '", xeniumOut, "' directory was not found.",
                "\n Verify 'xeniumOut' input directory.",
                call. = FALSE
            )
        resources <- .file_for_format(xeniumOut, format, ...)
    } else {
        stopifnot(
            (isScalarCharacter(resources) && file.exists(resources)) ||
                is(resources, "TENxFileList_OR_TENxH5")
        )
        if (
            !is(resources, "TENxFileList_OR_TENxH5") &&
            identical(tools::file_ext(resources), "h5")
        )
            resources <- TENxH5(resources, ...)
        else if (isScalarCharacter(resources))
            resources <- TENxFileList(resources, ...)
        xeniumOut <- dirname(path(resources))
    }

    xeniumfile <- .filter_xenium_file(xeniumOut)
    parqfile <- .filter_parquet_file(xeniumOut)
    bounds <- .boundaries_for_format(xeniumOut, boundaries_format)

    .TENxXenium(
        resources = resources,
        boundaries = bounds,
        coordNames = spatialCoordsNames,
        sampleId = sample_id,
        colData = parqfile,
        metadata = xeniumfile
    )
}

.validTENxXenium <- function(object) {
    isFL <- is(object@resources, "TENxFileList_OR_TENxH5")
    isSP <- is(object@colData, "TENxSpatialParquet")
    isXF <- is(object@metadata, "XeniumFile")
    if (all(isFL, isSP, isXF))
        TRUE
    else if (!isFL)
        "'resources' component is not of TENxFileList or TENxH5 class"
    else if (!isSP)
        "'colData' component is not of TENxSpatialParquet class"
    else if (!isXF)
        "'metadata' component is not of XeniumFile class"
}

#' @importFrom S4Vectors setValidity2
setValidity2("TENxXenium", .validTENxXenium)

# import TENxXenium method ------------------------------------------------

#' @describeIn TENxXenium-class Import Xenium Analyzer data
#'
#' @inheritParams BiocIO::import
#'
#' @importFrom BiocIO import
#' @importFrom methods as
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SummarizedExperiment assay assays rowData colData
#' @importFrom SingleCellExperiment mainExpName altExps
#'
#' @exportMethod import
setMethod("import", "TENxXenium", function(con, format, text, ...) {
    sce <- import(con@resources, ...)
    metadata <- import(con@metadata)
    coldata <- import(con@colData)

    ## TODO: mainExpName and altExpNames are lost when SCE sent to constructor
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = assay(sce)),
        rowData = rowData(sce),
        mainExpName = mainExpName(sce),
        altExps = altExps(sce),
        sample_id = con@sampleId,
        colData = as(coldata, "DataFrame"),
        spatialCoordsNames = con@coordNames,
        metadata = list(
            experiment.xenium = metadata,
            polygons = import(con@boundaries)
        )
    )
})
