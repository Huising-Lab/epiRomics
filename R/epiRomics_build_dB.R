#' Build epiRomics database
#'
#' @param epiRomics_db_file character string of path to properly formatted csv file containing epigenetic data. [See vignette for more details]
#' @param txdb_organism a character string containing the TxDB associated with your data.
#' @param epiRomics_genome a character string containing the genome associated with your data. e.g. "mm10" or "hg19".
#' @param epiRomics_organism a character string containing the org.db associated with your data.
#' @return Variable of class epiRomics for further downstream analysis
#' @export
#' @importFrom methods new
#' @importFrom data.table %like%
epiRomics_build_dB <-
  function(epiRomics_db_file,
           txdb_organism,
           epiRomics_genome,
           epiRomics_organism) {
    epiRomics_db_data_table <- utils::read.csv(epiRomics_db_file)

    for (i in 1:base::dim(epiRomics_db_data_table)[[1]]) {
      base::eval(base::parse(
        text = base::paste0(
          epiRomics_db_data_table[i, "name"],
          "_file <- ",
          "'",
          epiRomics_db_data_table[i, "path"] ,
          "'"
        )
      ))
      annotatr::read_annotations(
        con = base::eval(base::parse(text = (
          base::paste0(epiRomics_db_data_table[i, "name"],
                       "_file")
        ))),
        genome = epiRomics_db_data_table[i, "genome"],
        name = epiRomics_db_data_table[i, "name"],
        format = epiRomics_db_data_table[i, "format"]
      )
    }
    epiRomics_db_annot_list <-  base::c(# Pre-built (needs data.table package exported)
      annotatr::builtin_annotations()[annotatr::builtin_annotations() %like%
                                        epiRomics_genome],
      (annotatr::annotatr_cache$list_env()))
    epiRomics_dB <- methods::new("epiRomicsS4")

    epiRomics_dB@annotations <-
      annotatr::build_annotations(genome = epiRomics_genome, annotations = epiRomics_db_annot_list)
    # Fix later, too slow
    #    epiRomics_dB_meta <-
    #      base::data.frame(base::matrix("genomic", base::dim(base::as.data.frame(epiRomics_dB))[[1]], 3))
    #    base::colnames(epiRomics_dB_meta) <- base::c("name", "source", "type")
    #    for (i in 1:base::dim(base::as.data.frame(epiRomics_dB))[[1]]) {
    #      epiRomics_dB_meta[i,"name"] <- base::unlist(base::strsplit(base::as.data.frame(epiRomics_dB)[i,"id"], ":"))[[1]]
    #    }
    epiRomics_dB@meta <- epiRomics_db_data_table
    epiRomics_dB@txdb <- txdb_organism
    epiRomics_dB@organism <- epiRomics_organism
    epiRomics_dB@genome <-  epiRomics_genome
    base::return(epiRomics_dB)
  }
