# Public accessors for epiRomicsS4
#
# Getter and setter methods for every slot of the epiRomicsS4 class so
# users never need to reach into @slot or methods::slot().
#
# organism() extends the generic exported by BiocGenerics and genome()
# extends the generic exported by GenomeInfoDb, following Bioconductor
# convention. annotations(), meta(), and txdb() are new generics scoped
# to this package.

#' @include class-epiRomics.R
#' @importFrom BiocGenerics organism
#' @importFrom BiocGenerics "organism<-"
#' @importFrom GenomeInfoDb genome
#' @importFrom GenomeInfoDb "genome<-"
NULL

#' Accessors for the epiRomicsS4 class
#'
#' The \code{epiRomicsS4} class has five slots. Every slot is exposed
#' through a public getter and setter method so users should never need
#' to reach into the object with \code{obj@slot} or
#' \code{methods::slot(obj, "slot")}.
#'
#' \tabular{lll}{
#'   Slot \tab Getter \tab Setter \cr
#'   \code{annotations} \tab \code{\link{annotations}(x)} \tab
#'     \code{\link[=annotations<-]{annotations(x) <- value}} \cr
#'   \code{meta} \tab \code{\link{meta}(x)} \tab
#'     \code{\link[=meta<-]{meta(x) <- value}} \cr
#'   \code{txdb} \tab \code{\link{txdb}(x)} \tab
#'     \code{\link[=txdb<-]{txdb(x) <- value}} \cr
#'   \code{organism} \tab \code{\link{organism}(object)} \tab
#'     \code{\link[=organism<-]{organism(object) <- value}} \cr
#'   \code{genome} \tab \code{\link{genome}(x)} \tab
#'     \code{\link[=genome<-]{genome(x) <- value}} \cr
#' }
#'
#' \code{organism()} extends the generic from \pkg{BiocGenerics} and
#' \code{genome()} extends the generic from \pkg{GenomeInfoDb}, so
#' \code{epiRomicsS4} objects respond to those accessors the same way
#' other Bioconductor objects do. \code{annotations()}, \code{meta()},
#' and \code{txdb()} are generics scoped to \pkg{epiRomics}.
#'
#' Every setter validates the updated object via
#' \code{\link[methods]{validObject}} before returning it, so invalid
#' assignments (for example an empty-string \code{genome}) fail fast
#' with an informative error.
#'
#' @return An overview topic; see the individual accessor pages for
#'   call signatures and return values.
#' @examples
#' db <- make_example_database()
#'
#' # Read every slot through its getter
#' annotations(db)
#' meta(db)
#' txdb(db)
#' organism(db)
#' genome(db)
#'
#' # Setters return an updated object
#' db2 <- db
#' genome(db2) <- "mm10"
#' organism(db2) <- "org.Mm.eg.db"
#' genome(db2)
#' organism(db2)
#' @seealso \code{\linkS4class{epiRomicsS4}}
#' @name epiRomicsS4-accessors
#' @aliases epiRomicsS4-accessors
#' @docType methods
NULL

# ---------------------------------------------------------------------------
# annotations() getter
# ---------------------------------------------------------------------------

#' Access the genomic annotations slot of an epiRomicsS4 object
#'
#' Returns the \code{GRanges} object stored in the \code{annotations}
#' slot. Public accessor replacement for \code{obj@annotations} or
#' \code{methods::slot(obj, "annotations")}.
#'
#' @param x An \code{\linkS4class{epiRomicsS4}} object.
#' @param ... Currently unused; reserved for method extension.
#' @return A \code{\link[GenomicRanges]{GRanges}} object. Empty
#'   (length 0) if no annotations have been assigned.
#' @examples
#' db <- make_example_database()
#' annotations(db)
#' length(annotations(db))
#' @seealso \code{\link{annotations<-}}, \code{\linkS4class{epiRomicsS4}}
#' @export
methods::setGeneric(
  "annotations",
  function(x, ...) standardGeneric("annotations")
)

#' @rdname annotations
#' @export
methods::setMethod(
  "annotations", "epiRomicsS4",
  function(x, ...) x@annotations
)

# ---------------------------------------------------------------------------
# annotations<-() setter
# ---------------------------------------------------------------------------

#' Replace the annotations slot of an epiRomicsS4 object
#'
#' Assigns a new \code{GRanges} object to the \code{annotations} slot.
#' Validity is checked via \code{\link[methods]{validObject}}.
#'
#' @param x An \code{\linkS4class{epiRomicsS4}} object.
#' @param ... Currently unused; reserved for method extension.
#' @param value A \code{\link[GenomicRanges]{GRanges}} object.
#' @return The updated \code{\linkS4class{epiRomicsS4}} object.
#' @examples
#' db <- make_example_database()
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
#' gr$type <- "hg38_custom_h3k4me1"
#' annotations(db) <- gr
#' length(annotations(db))
#' @name annotations-set
#' @aliases annotations<- annotations<-,epiRomicsS4-method
#' @export
methods::setGeneric(
  "annotations<-",
  function(x, ..., value) standardGeneric("annotations<-")
)

#' @rdname annotations-set
#' @export
methods::setMethod(
  "annotations<-", "epiRomicsS4",
  function(x, ..., value) {
    x@annotations <- value
    methods::validObject(x)
    x
  }
)

# ---------------------------------------------------------------------------
# meta() getter
# ---------------------------------------------------------------------------

#' Access the metadata slot of an epiRomicsS4 object
#'
#' Returns the \code{data.frame} stored in the \code{meta} slot, listing
#' the name and type of every data source loaded into the database.
#' Public accessor replacement for \code{obj@meta} or
#' \code{methods::slot(obj, "meta")}.
#'
#' @param x An \code{\linkS4class{epiRomicsS4}} object.
#' @param ... Currently unused; reserved for method extension.
#' @return A \code{data.frame} with one row per data source. Empty
#'   (zero rows) if no metadata has been assigned.
#' @examples
#' db <- make_example_database()
#' meta(db)
#' nrow(meta(db))
#' @seealso \code{\link{meta<-}}, \code{\linkS4class{epiRomicsS4}}
#' @export
methods::setGeneric(
  "meta",
  function(x, ...) standardGeneric("meta")
)

#' @rdname meta
#' @export
methods::setMethod(
  "meta", "epiRomicsS4",
  function(x, ...) x@meta
)

# ---------------------------------------------------------------------------
# meta<-() setter
# ---------------------------------------------------------------------------

#' Replace the metadata slot of an epiRomicsS4 object
#'
#' Assigns a new \code{data.frame} to the \code{meta} slot.
#'
#' @param x An \code{\linkS4class{epiRomicsS4}} object.
#' @param ... Currently unused; reserved for method extension.
#' @param value A \code{data.frame}.
#' @return The updated \code{\linkS4class{epiRomicsS4}} object.
#' @examples
#' db <- make_example_database()
#' new_meta <- meta(db)
#' new_meta$extra <- "annotated"
#' meta(db) <- new_meta
#' colnames(meta(db))
#' @name meta-set
#' @aliases meta<- meta<-,epiRomicsS4-method
#' @export
methods::setGeneric(
  "meta<-",
  function(x, ..., value) standardGeneric("meta<-")
)

#' @rdname meta-set
#' @export
methods::setMethod(
  "meta<-", "epiRomicsS4",
  function(x, ..., value) {
    x@meta <- value
    methods::validObject(x)
    x
  }
)

# ---------------------------------------------------------------------------
# txdb() getter
# ---------------------------------------------------------------------------

#' Access the TxDb package::object identifier slot
#'
#' Returns the character string stored in the \code{txdb} slot. The value
#' has the form \code{"TxDb.Package::TxDb.object"}, e.g.
#' \code{"TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene"}.
#' Public accessor replacement for \code{obj@txdb} or
#' \code{methods::slot(obj, "txdb")}.
#'
#' @param x An \code{\linkS4class{epiRomicsS4}} object.
#' @param ... Currently unused; reserved for method extension.
#' @return A character scalar. Empty character if unset.
#' @examples
#' db <- make_example_database()
#' txdb(db)
#' @seealso \code{\link{txdb<-}}, \code{\linkS4class{epiRomicsS4}}
#' @export
methods::setGeneric(
  "txdb",
  function(x, ...) standardGeneric("txdb")
)

#' @rdname txdb
#' @export
methods::setMethod(
  "txdb", "epiRomicsS4",
  function(x, ...) x@txdb
)

# ---------------------------------------------------------------------------
# txdb<-() setter
# ---------------------------------------------------------------------------

#' Replace the TxDb identifier slot of an epiRomicsS4 object
#'
#' Assigns a new \code{"TxDb.Package::TxDb.object"} character string to
#' the \code{txdb} slot. Validity is checked via
#' \code{\link[methods]{validObject}}.
#'
#' @param x An \code{\linkS4class{epiRomicsS4}} object.
#' @param ... Currently unused; reserved for method extension.
#' @param value A character scalar in
#'   \code{"TxDb.Package::TxDb.object"} form.
#' @return The updated \code{\linkS4class{epiRomicsS4}} object.
#' @examples
#' db <- make_example_database()
#' txdb(db) <- "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene"
#' txdb(db)
#' @name txdb-set
#' @aliases txdb<- txdb<-,epiRomicsS4-method
#' @export
methods::setGeneric(
  "txdb<-",
  function(x, ..., value) standardGeneric("txdb<-")
)

#' @rdname txdb-set
#' @export
methods::setMethod(
  "txdb<-", "epiRomicsS4",
  function(x, ..., value) {
    x@txdb <- value
    methods::validObject(x)
    x
  }
)

# ---------------------------------------------------------------------------
# organism() method on the BiocGenerics generic
# ---------------------------------------------------------------------------

#' Access the organism annotation package name
#'
#' Returns the character string stored in the \code{organism} slot
#' (e.g. \code{"org.Hs.eg.db"}, \code{"org.Mm.eg.db"}). This method
#' extends the \code{organism()} generic from the
#' \pkg{BiocGenerics} package.
#'
#' @param object An \code{\linkS4class{epiRomicsS4}} object.
#' @return A character scalar. Empty character if unset.
#' @examples
#' db <- make_example_database()
#' organism(db)
#' @name organism
#' @aliases organism,epiRomicsS4-method
#' @exportMethod organism
methods::setMethod(
  "organism", "epiRomicsS4",
  function(object) object@organism
)

#' Replace the organism annotation package name
#'
#' Assigns a new \code{org.*.eg.db} package name to the \code{organism}
#' slot. This method extends the \code{organism<-()} replacement
#' generic from the \pkg{BiocGenerics} package.
#'
#' @param object An \code{\linkS4class{epiRomicsS4}} object.
#' @param value A character scalar, typically an
#'   \code{org.*.eg.db} package name.
#' @return The updated \code{\linkS4class{epiRomicsS4}} object.
#' @examples
#' db <- make_example_database()
#' organism(db) <- "org.Hs.eg.db"
#' organism(db)
#' @name organism-set
#' @aliases organism<- organism<-,epiRomicsS4-method
#' @exportMethod "organism<-"
methods::setMethod(
  "organism<-", "epiRomicsS4",
  function(object, value) {
    object@organism <- value
    methods::validObject(object)
    object
  }
)

# ---------------------------------------------------------------------------
# genome() method on the GenomeInfoDb generic
# ---------------------------------------------------------------------------

#' Access the genome assembly name
#'
#' Returns the character string stored in the \code{genome} slot
#' (e.g. \code{"hg38"}, \code{"mm10"}). This method extends the
#' \code{genome()} generic from the \pkg{GenomeInfoDb} package.
#'
#' @param x An \code{\linkS4class{epiRomicsS4}} object.
#' @return A character scalar. Empty character if unset.
#' @examples
#' db <- make_example_database()
#' genome(db)
#' @name genome
#' @aliases genome,epiRomicsS4-method
#' @exportMethod genome
methods::setMethod(
  "genome", "epiRomicsS4",
  function(x) x@genome
)

#' Replace the genome assembly name
#'
#' Assigns a new genome assembly name (e.g. \code{"hg38"}, \code{"mm10"})
#' to the \code{genome} slot. This method extends the
#' \code{genome<-()} replacement generic from the \pkg{GenomeInfoDb}
#' package.
#'
#' @param x An \code{\linkS4class{epiRomicsS4}} object.
#' @param value A non-empty character scalar.
#' @return The updated \code{\linkS4class{epiRomicsS4}} object.
#' @examples
#' db <- make_example_database()
#' genome(db) <- "hg38"
#' genome(db)
#' @name genome-set
#' @aliases genome<- genome<-,epiRomicsS4-method
#' @exportMethod "genome<-"
methods::setMethod(
  "genome<-", "epiRomicsS4",
  function(x, value) {
    x@genome <- value
    methods::validObject(x)
    x
  }
)
