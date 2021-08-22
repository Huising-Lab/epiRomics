#' zzz
#'
#' zzzz
#'
#' @param libname
#'
#' @param pkgname
#'
#' @noRd

.onLoad <-
  function(libname = find.package("epiRomics"),
           pkgname = "epiRomics") {
    # do whatever needs to be done when the package is loaded
    # some people use it to bombard users with
    # messages using
    base::packageStartupMessage("epiRomics package loaded.")
    base::packageStartupMessage("Automatic dependency check and verification of sample data presence")

    options(timeout = max(3000, getOption("timeout")))

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    BiocManager::install()


    to_install_cran <-
      c(
        "data.table",
        "party",
        "plyr"
      )
    for (i in to_install_cran) {
      base::packageStartupMessage(paste("looking for ", i))
      if (!requireNamespace(i)) {
        base::packageStartupMessage(paste("     installing", i))
        utils::install.packages(i)
      }
    }

    to_install_bc <-
      c(
        "AnnotationDbi",
        "annotatr",
        "BiocGenerics",
        "GenomicFeatures",
        "GenomicRanges",
        "Gviz",
        "IRanges",
        "rtracklayer",
        "org.Hs.eg.db",
        "TxDb.Hsapiens.UCSC.hg38.knownGene"
      )

    for (i in to_install_bc) {
      base::packageStartupMessage(paste("looking for ", i))
      if (!requireNamespace(i)) {
        base::packageStartupMessage(paste("     installing", i))
        BiocManager::install(i, type = "source")
      }
    }


    # Download ChIP, Histone, Functional, and DBA Alpha versus Beta

    # BED
    destfile <- paste0(
      system.file(package = "epiRomics"),
      "/extdata/BED_Annotation.tar.gz"
    )

    fileURL <-
      "https://dl.dropboxusercontent.com/s/qf0hca9rnrbn70j/BED_Annotation.tar.gz"
    if (!file.exists(destfile)) {
      utils::download.file(fileURL, destfile, method = "auto")
    }

    utils::untar(
      paste0(
        system.file(package = "epiRomics"),
        "/extdata/BED_Annotation.tar.gz"
      ),
      exdir = paste0(system.file(package = "epiRomics"), "/extdata/")
    )


    # CHIP
    destfile <- paste0(
      system.file(package = "epiRomics"),
      "/extdata/ChIP.tar.gz"
    )
    fileURL <-
      "https://dl.dropboxusercontent.com/s/uiv4wu45rvty8pt/ChIP.tar.gz"
    if (!file.exists(destfile)) {
      utils::download.file(fileURL, destfile, method = "auto")
    }

    utils::untar(paste0(system.file(package = "epiRomics"), "/extdata/ChIP.tar.gz"),
      exdir = paste0(system.file(package = "epiRomics"), "/extdata/")
    )



    # Histone
    destfile <- paste0(
      system.file(package = "epiRomics"),
      "/extdata/Histone.tar.gz"
    )
    fileURL <-
      "https://dl.dropboxusercontent.com/s/7zf6ufhswwe6so5/Histone.tar.gz"
    if (!file.exists(destfile)) {
      utils::download.file(fileURL, destfile, method = "auto")
    }

    utils::untar(
      paste0(
        system.file(package = "epiRomics"),
        "/extdata/Histone.tar.gz"
      ),
      exdir = paste0(system.file(package = "epiRomics"), "/extdata/")
    )

    # Bigwig
    destfile <- paste0(
      system.file(package = "epiRomics"),
      "/extdata/BigWigs.tar.gz"
    )
    fileURL <-
      "https://dl.dropboxusercontent.com/s/qwp3pt5d9ot60nr/BigWigs.tar.gz"
    if (!file.exists(destfile)) {
      utils::download.file(fileURL, destfile, method = "auto")
    }
    utils::untar(
      paste0(
        system.file(package = "epiRomics"),
        "/extdata/BigWigs.tar.gz"
      ),
      exdir = paste0(system.file(package = "epiRomics"), "/extdata/")
    )

    # DBA
    destfile <- paste0(
      system.file(package = "epiRomics"),
      "/extdata/DBA_Beta_Versus_Alpha.tar.gz"
    )
    fileURL <-
      "https://dl.dropboxusercontent.com/s/28pyns5ktzzggmm/DBA_Beta_Versus_Alpha.tar.gz"
    if (!file.exists(destfile)) {
      utils::download.file(fileURL, destfile, method = "auto")
    }
    utils::untar(
      paste0(
        system.file(package = "epiRomics"),
        "/extdata/DBA_Beta_Versus_Alpha.tar.gz"
      ),
      exdir = paste0(system.file(package = "epiRomics"), "/extdata/")
    )



    # Fix paths for dB build
    example_epiRomics_Db_sheet <-
      utils::read.csv(file = system.file("extdata", "example_epiRomics_Db_sheet.csv", package = "epiRomics"))
    example_epiRomics_Db_sheet$path <-
      paste0(
        system.file(package = "epiRomics"),
        "/extdata/",
        example_epiRomics_Db_sheet$path
      )
    data.table::fwrite(
      example_epiRomics_Db_sheet,
      file = paste0(
        system.file(package = "epiRomics"),
        "/extdata/example_epiRomics_Db_sheet_user_paths.csv"
      ),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )




    # Fix paths for bigwig connection
    epiRomics_track_connection <-
      utils::read.csv(system.file("extdata", "example_epiRomics_BW_sheet.csv", package = "epiRomics"))
    epiRomics_track_connection$path <-
      paste0(
        system.file(package = "epiRomics"),
        "/extdata/",
        epiRomics_track_connection$path
      )


    data.table::fwrite(
      epiRomics_track_connection,
      file = paste0(
        system.file(package = "epiRomics"),
        "/extdata/example_epiRomics_BW_sheet_user_paths.csv"
      ),
      sep = ",",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )




    base::packageStartupMessage("You are ready to go. For feedback, please email: ammawla@ucdavis.edu")
  }
