#' Predicts TF behavior in association with enhanceosome presence in accordance with functional annotation
#'
#' @param epiRomics_putative_enhanceosome epiRomics class database containing putative enhanceosome calls
#' @return Returned decision tree available for plotting
#' @export


epiRomics_predictors <- function(epiRomics_putative_enhanceosome) {
  epiRomics_data_test <- base::as.data.frame(epiRomics_putative_enhanceosome@annotations)
  epiRomics_chip_names <- epiRomics_putative_enhanceosome@meta$name[epiRomics_putative_enhanceosome@meta$type ==
    "chip"]
  epiRomics_combos <- base::as.data.frame(matrix(nrow = dim(epiRomics_data_test)[1], ncol = 0))
  for (i in 1:(base::length(epiRomics_chip_names) - 1)) {
    n <- i + 1
    for (j in n:(base::length(epiRomics_chip_names) - 1)) {
      epiRomics_combos <- base::cbind(epiRomics_combos, base::eval(base::parse(text = "epiRomics_data_test[, epiRomics_chip_names[i]] * epiRomics_data_test[, epiRomics_chip_names[j]]")))
      base::colnames(epiRomics_combos)[base::dim(epiRomics_combos)[2]] <- base::paste0(epiRomics_chip_names[i],
        sep = "_", epiRomics_chip_names[j]
      )
    }
  }

  # epiRomics_data_test <- cbind(epiRomics_combos,epiRomics_data_test[, c(epiRomics_chip_names,
  # 'annotation')])
  epiRomics_data_test <- epiRomics_data_test[, c(epiRomics_chip_names, "annotation")]
  epiRomics_annotation_filter <- epiRomics_data_test$annotation
  for (i in 1:base::length(epiRomics_annotation_filter)) {
    epiRomics_annotation_filter[i] <- base::unlist(base::strsplit(
      x = epiRomics_annotation_filter[i],
      split = "(", fixed = TRUE
    ))[[1]]
  }
  epiRomics_data_test$annotation <- epiRomics_annotation_filter
  epiRomics_data_test$annotation[epiRomics_data_test$annotation != "Distal Intergenic"] <- "Genic"
  epiRomics_data_test$annotation <- base::as.factor(epiRomics_data_test$annotation)
  output.tree <- party::ctree(annotation ~ ., data = epiRomics_data_test)
  # Number of co-TFs binding determines functional genome table(predict(output.tree),
  # epiRomics_data_test$ChIP_Hits)
  base::return(output.tree)
}
