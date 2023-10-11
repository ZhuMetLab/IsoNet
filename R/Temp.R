################################################################################
# setting style ----------------------------------------------------------------
errorStyle <- crayon::red
warnStyle <- crayon::magenta
spliteStyle <- crayon::silver
startStyle <- crayon::blue
stageStyle <- crayon::yellow
normalStyle <- crayon::green


################################################################################
# 20220901 edge rule -----------------------------------------------------------
# mass_diff_table <- readr::read_csv('G:/00_projects/06_isotopeTracing/01_package/data/20220728_edge_rules.csv')
#
# idx <- which(mass_diff_table$AtomicDifference %in% c('none_adduct',
#                                                      'transformation',
#                                                      'unspecified',
#                                                      'unspecified|unspecified|unspecified',
#                                                      'unspecified|unspecified|unspecified|unspecified|unspecified'))
# mass_diff_table <- mass_diff_table[-idx, ]
#
# mass_diff_table <- mass_diff_table %>%
#   dplyr::mutate(GroupDeltaMass_abs = abs(GroupDeltaMass))
#
#
# usethis::use_data(mass_diff_table, overwrite = TRUE)
