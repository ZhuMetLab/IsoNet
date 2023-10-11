# constructMIDNet(tracer_table_file = 'MID_more_than_3_POS_maxC.csv',
#                 ms2_file = 'ms2_data.msp',
#                 sample_names = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6"),
#                 dir_path = "G:/00_projects/06_isotopeTracing/data_processing/20220906_MIDNet",
#                 mid_cutoff = 0.95,
#                 mid_fc = 20,
#                 mid_isoDegree = 0.02,
#                 mid_min_motifLen = "50%",
#                 mid_max_motif = 1000,
#                 scoring_approach = 'gnps',
#                 mz_tol_ms2 = 15,
#                 ms2_score_cutoff = 0.5,
#                 num_matched_frag_cutoff = 1)
#
#
################################################################################
#
# load('G:/00_projects/06_isotopeTracing/data_processing/20220831_MIDNet/00_intermidate_data/list_SIMres.RData')
#
# list_SIMres
# which(NodeTable$From == 'M292T447NA5' & NodeTable$To == 'M380T404C03451')
#
# NodeTable[501567,]
# NodeTable[248382,]

# library(CalMIDSIM)
# setwd("D:/gaoyang/MID_caculation_test/20220719_POS_293T")
#
# calMIDSim(NodeTable = NULL,
#           TargetTable_dir = "G:/00_projects/06_isotopeTracing/test/MID_more_than_3_POS_maxC.csv",
#           TargetSample = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6"),
#           Cutoff = 0.95,
#           FoldChange = 20,
#           IsoDegree = 0.02,
#           MinMotifLength = "60%",
#           printall = FALSE,
#           maxMotif = 1000)
#
#
# NodeTable = NULL
# Pattern_name = "Unknown"
# TargetTable_dir = 'G:/00_projects/06_isotopeTracing/data_processing/20220831_MIDNet'
# TargetSample = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6")
# Cutoff = 0.9
# IsoDegree = 0.1
# FoldChange = 20
# MinMotifLength = "50%"
# printall = FALSE
# maxMotif = 10000
#
#
#
# MIDSIMCal(TargetTable_dir = "MID_more_than_3_POS_maxC.csv",
#           TargetSample = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6"),
#           Cutoff = 0.95,
#           FoldChange = 20,
#           IsoDegree = 0.02,
#           MinMotifLength = "50%",
#           printall = FALSE,
#           maxMotif = 1000)
#
#
# constructMIDNet(tracer_table_file = "MID_more_than_3_POS_maxC.csv",
#                 ms2_file = "ms2_data.msp",
#                 sample_names = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6"),
#                 dir_path = 'G:/00_projects/06_isotopeTracing/data_processing/20220831_MIDNet',
#                 mid_cutoff = 0.95,
#                 mid_fc = 20,
#                 mid_isoDegree = 0.02,
#                 mid_min_motifLen = "50%",
#                 mid_max_motif = 1000)
#
#
# test <- c(pairs$From_feature, pairs$To_feature) %>% unique()
# temp_idx <- which(!(test %in% names(ms2_spec)))
# test <- test[temp_idx]
# pairs <- pairs %>%
#   dplyr::filter(!(From_feature %in% test | To_feature %in% test))
#
#
#
################################################################################
#
# constructMIDNet(tracer_table_file = "MID_more_than_3_POS_maxC.csv",
#                 ms2_file = "ms2_data.msp",
#                 sample_names = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6"),
#                 dir_path = 'G:/00_projects/06_isotopeTracing/data_processing/20220831_MIDNet',
#                 multi_thread = 6,
#                 mid_cutoff = 0.95,
#                 mid_fc = 20,
#                 mid_isoDegree = 0.02,
#                 mid_min_motifLen = "50%",
#                 mid_max_motif = 1000)
#
# constructMIDNet(tracer_table_file = "MID_more_than_3_POS_maxC.csv",
#                 ms2_file = "ms2_data.msp",
#                 sample_names = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6"),
#                 dir_path = 'G:/00_projects/06_isotopeTracing/data_processing/20220923_MIDNet/',
#                 multi_thread = 6,
#                 mid_cutoff = 0.95,
#                 mid_fc = 20,
#                 mid_isoDegree = 0.02,
#                 mid_min_motifLen = "50%",
#                 mid_max_motif = 1000)
