#!/usr/bin/env Rscript
parse_args = function(x) {
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)
  # giving errors when we have lists like c(xxx, xxx) since it will separate it
  # args_list = unlist(strsplit(x, ', ')[[1]])
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
  # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
  args_vals = lapply(args_list, function(x) {
    x_splt = strsplit(x, split=":")[[1]]
    c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
  })

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

  parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}
opt = list(
    prefix = ifelse("$task.ext.prefix" == "null", "$meta.id", "$task.ext.prefix")
)
args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]
print(opt)


# Script ####
library(dplyr)
library(CNAqc)

multi_cnaqc = readRDS(file = "$rds_join")
if (as.logical(opt[["qc_chr"]])){
  multisample_jointTable <- lapply(names(multi_cnaqc), function(c){
    mutations_multisample_chr = get_sample(m_cnaqc_obj = multi_cnaqc[[c]],sample = get_sample_name(multi_cnaqc[[c]]),
                                           which_obj = "original")
    multisample_jointTable_chr = list()
    for (s in get_sample_name(multi_cnaqc[[c]])){
      purity = mutations_multisample_chr[[s]][["purity"]]
      multisample_jointTable_chr[[s]] = mutations_multisample_chr[[s]][["mutations"]] %>%
        dplyr::mutate(purity = purity) %>%
        dplyr::mutate(patient_id = "$meta.patient")
    }
    joint_table = bind_rows(multisample_jointTable_chr)
    return(joint_table)
  })
} else{
  mutations_multisample <- get_sample(m_cnaqc_obj = multi_cnaqc,sample = get_sample_name(multi_cnaqc),
                                      which_obj = "original")
  multisample_jointTable = list()

  for (s in get_sample_name(multi_cnaqc)){
    purity = mutations_multisample[[s]][["purity"]]
    multisample_jointTable[[s]] = mutations_multisample[[s]][["mutations"]] %>%
      dplyr::mutate(purity = purity) %>%
      dplyr::mutate(patient_id = "$meta.patient")

    if ('CCF_estimates' %in%  names(mutations_multisample[[s]])){

      ccf = lapply(mutations_multisample[[s]][["CCF_estimates"]], FUN = function(k){
        k[["mutations"]] %>% select(mutation_id, CCF, mutation_multiplicity)
      }) %>% bind_rows()

      multisample_jointTable[[s]] = multisample_jointTable[[s]] %>% left_join(ccf, by = join_by(mutation_id))
    }
  }
}


joint_table = bind_rows(multisample_jointTable)

joint_table <- joint_table %>%
  dplyr::mutate(
    from = as.integer(from),
    to   = as.integer(from + nchar(ref) - 1)
  )

options(scipen = 999)

write.table(joint_table, file = paste0(opt[["prefix"]], "_joint_table.tsv"), append = F, quote = F, sep = "\t", row.names = FALSE)

# version export
f = file("versions.yml","w")
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
close(f)
