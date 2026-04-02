#!/usr/bin/env Rscript

opt = list(
    prefix = ifelse("$task.ext.prefix" == "null", "$meta.id", "$task.ext.prefix")
)

# Script ####
library(dplyr)
library(CNAqc)

multi_cnaqc = readRDS(file = "$rds_join")
mutations_multisample = get_sample(m_cnaqc_obj = multi_cnaqc,sample = get_sample_name(multi_cnaqc),
                                   which_obj = "original")
multisample_jointTable = list()

for (s in get_sample_name(multi_cnaqc)){
    purity = mutations_multisample[[s]][["purity"]]
    multisample_jointTable[[s]] = mutations_multisample[[s]][["mutations"]] %>%
      dplyr::mutate(purity = purity) %>%
      dplyr::mutate(patient_id = "$meta.patient")
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
