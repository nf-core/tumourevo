get_sample = function(m_cnaqc_obj,
                      sample,
                      which_obj) {
  if (class(m_cnaqc_obj) != "m_cnaqc") {
    cli::cli_abort(c(
      "mCNAqc object was build without keeping original CNAqc objects",
      "x" = "It is not possible to retrieve the required samples"
    ))
  }

  consented_obj = c("shared", "original")

  if ((which_obj %in% consented_obj) == FALSE) {
    cli::cli_abort("{.var which_obj} must be one of 'shared' or 'original'")
  }
  # define the element names
  if (which_obj == "original") {
    type = "original_cnaqc_objc"
    # check if the original cnaqc obj exist
    check_or = any(names(m_cnaqc_obj) == type)
    if (check_or == FALSE) {
      cli::cli_abort(c(
        "mCNAqc object was build without keeping original CNAqc objects",
        "x" = "It is not possible to retrieve the required samples"
      ))
    } else {
      cli::cli_h1("Retrieving original {.cls CNAqc} objects")
      cnaqc_samples = m_cnaqc_obj[[type]][sample]
    }
  } else {
    type = "cnaqc_obj_new_segmentation"
    cli::cli_h1("Retrieving {.cls CNAqc} objects with the new segmentation")
    cnaqc_samples = m_cnaqc_obj[[type]][sample]
  }
  return(cnaqc_samples)
}

get_sample_name = function(x) {
  if (class(x) == "m_cnaqc") {
    lapply(x$cnaqc_obj_new_segmentation, function(y) {
      y$sample
    }) %>% unlist() %>% unname()

  } else if (class(x) == "cnaqc") {
    x$sample
  } else {
    cli::cli_abort(c(
      "Must provide a {.field m_cnaqc} object",
      "x" = "{.var x} is a {.cls {class(x)}}"
    ))
  }
}

get_mCNAqc_stats = function(m_cnaqc_obj) {
  stats = m_cnaqc_obj[["m_cnaqc_stats"]]
  return(stats)
}
