################################################################################
#' @title constructMIDNet
#' @author Zhiwei Zhou
#' @param tracer_table_file
#' @param ms2_file
#' @param dir_path
#' @param mid_cutoff
#' @param mid_fc
#' @param mid_isoDegree
#' @param mid_min_motifLen
#' @param min_max_motif
#' @param scoring_approach
#' @param mz_tol_ms2
#' @param mz_tol_label
#' @param mass_diff_freq_cut_off
#' @import crayon
#' @importFrom magrittr '%>%'
#' @export

# c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6")

constructMIDNet <- function(tracer_table_file = "MID_more_than_3_POS_maxC.csv",
                            ms2_file = "ms2_data.msp",
                            id_table = 'id_table.csv',
                            sample_names,
                            dir_path = '.',
                            multi_thread = 1,

                            # parameters: calculating MID similarity
                            mid_cutoff = 0.95,
                            mid_fc = 20,
                            mid_isoDegree = 0.02,
                            mid_min_motifLen = 0.6,
                            mid_max_motif = 1000,
                            ignore_max = FALSE,

                            # parameters: calculating ms2 network
                            scoring_approach = 'gnps',
                            mz_tol_ms2 = 15,
                            ms2_score_cutoff = 0.5,
                            num_matched_frag_cutoff = 1,

                            # parameters: add edge labels
                            mz_tol_label = 15,
                            mass_diff_freq_cut_off = 1
                            ) {
  cat(startStyle("IsoNet run log\n"))
  cat(startStyle(as.character(Sys.time()), '\n\n'))

  # create run.log file
  cat("IsoNet run log\n", file = file.path(dir_path, "run.log.txt"))
  cat('IsoNet verison:', as.character(packageVersion('IsoNet')), '\n', file = file.path(dir_path, "run.log.txt"), append = TRUE)
  cat(as.character(Sys.time()), '\n', file = file.path(dir_path, "run.log.txt"), append = TRUE)


  # check file existed
  cat(spliteStyle("===============================================\n"))
  cat(stageStyle("Step 1: check required files ...\n"))
  cat("Step 1: check required files ...\n", file = file.path(dir_path, "run.log.txt"), append = TRUE)
  checkFiles(dir_path = dir_path,
             tracer_table_file = tracer_table_file,
             ms2_file = ms2_file,
             id_table = id_table)

  # construct MID network
  cat(spliteStyle("===============================================\n"))
  cat(stageStyle("Step 2: calculate isotopologue similarity ...\n"))
  cat("Step 2: calculate isotopologue similarity ...\n", file = file.path(dir_path, "run.log.txt"), append = TRUE)
  if (!('GlobalMatchResult_filterred.csv' %in% list.files(file.path(dir_path, '00_intermidate_data')))) {
    calMIDSim(NodeTable = NULL,
              TargetTable_dir = file.path('00_intermidate_data', 'trace_table_clean.csv'),
              TargetSample = sample_names,
              dir_path = dir_path,
              multi_thread = multi_thread,
              Cutoff = mid_cutoff,
              FoldChange = mid_fc,
              IsoDegree = mid_isoDegree,
              MinMotifLength = mid_min_motifLen,
              maxMotif = mid_max_motif,
              ignore_max = ignore_max)
  } else {
    cat(warnStyle('Warning: the MID similarity file existed, skip this step!\n'))
  }


  # construct MS2 network
  cat(spliteStyle("===============================================\n"))
  cat(stageStyle("Step 3: calculate molecule network similarity ...\n"))
  cat("Step 3: calculate molecule network similarity ...\n", file = file.path(dir_path, "run.log.txt"), append = TRUE)

  if (!('ms2_sim_table.csv' %in% list.files(file.path(dir_path, '00_intermidate_data')))) {
    calMS2Net(network_file = 'GlobalMatchResult_filterred.csv',
              dir_path = dir_path,
              ms2_file = ms2_file,
              scoring_approach = scoring_approach,
              mz_tol_ms2 = mz_tol_ms2,
              ms2_score_cutoff = ms2_score_cutoff,
              num_matched_frag_cutoff = num_matched_frag_cutoff)
  } else {
    cat(warnStyle('Warning: the MS2 similarity file existed, skip this step!\n'))
  }


  # output
  cat(spliteStyle("===============================================\n"))
  cat(stageStyle("Step 4: export MID networks ...\n"))
  cat("Step 4: export MID networks ...\n", file = file.path(dir_path, "run.log.txt"), append = TRUE)

  # Node table
  addNodeLabel(tracer_table_file = 'trace_table_clean.csv',
               network_file = 'ms2_sim_table_filter.csv',
               dir_path = dir_path)

  # Edge table
  addEdgeLabel(tracer_table_file = 'trace_table_clean.csv',
               network_file = 'ms2_sim_table_filter.csv',
               id_table = id_table,
               dir_path = dir_path,
               mz_ppm = mz_tol_label,
               mass_diff_freq_cut_off = mass_diff_freq_cut_off)

  # Add ID and Step
  calculateStep(id_table = id_table,
                dir_path = dir_path)


  cat("===============================================\n")
  cat(stageStyle("Done!\n"))
  cat(startStyle(as.character(Sys.time()), '\n\n'))
  cat(as.character(Sys.time()), file = file.path(dir_path, "run.log.txt"), append = TRUE)

}

################################################################################
# checkFiles -------------------------------------------------------------------

checkFiles <- function(dir_path = '.',
                       tracer_table_file = "G:/00_projects/06_isotopeTracing/test/MID_more_than_3_POS_maxC.csv",
                       ms2_file = "ms2_data.msp",
                       id_table = id_table
){
  # browser()
  # check existences
  list_file_1 <- list.files(dir_path, recursive = F, full.names = F)
  list_file <- list.files(file.path(dir_path, 'Result'), recursive = F, full.names = T)
  if (!(tracer_table_file %in% list_file)) {
    stop('The tracer table', tracer_table_file, 'NOT existed!\n')
  }
  if (!(ms2_file %in% list_file_1)) {
    stop(errorStyle('The MS2 file', ms2_file, 'NOT existed!\n'))
  }

  if (!(id_table %in% list_file)) {
    stop(errorStyle('The id table ', id_table, 'NOT existed!\n'))
  }

  cat(normalStyle('Required files: Checked!\n'))

  # check format
    # tracer_table
  temp_data <- readr::read_csv(file.path(dir_path, tracer_table_file)
                               # , show_col_types = FALSE
                               )[,-1] # skip the first useless column (20230420)
  temp <- match(c('name', 'id', 'mz', 'formula', 'compound', 'label'), colnames(temp_data))
  if (!identical(temp, 1:6)) {
    stop(errorStyle('Plese check tracer_table_file format\n'))
  }

  if (!is.numeric(temp_data$label)) {
    stop(errorStyle('The label column is not numeric\n'))
  }
  cat(normalStyle('Tracer table format check: Checked!\n'))

  # ms2_file
  if (!(all(stringr::str_detect(ms2_file, pattern = '(\\.msp|\\.mgf)')))) {
    stop(errorStyle('IsoNet only support mgf/msp formats.'))
  }
  cat(normalStyle('MS2 file format check: Checked!\n'))

  # remove features without ms2
  ms2_spec <- readMSP(file.path(dir_path, ms2_file), mode = 'all')
  names(ms2_spec) <- sapply(ms2_spec, function(x){
    x$info$NAME
  })
  idx_rm <- which(!(temp_data$name %in% names(ms2_spec)))
  if (length(idx_rm) == 0) {
    cat(normalStyle('All features have MS/MS spectra: Checked!\n'))

    dir.create(file.path(dir_path, '00_intermidate_data'), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(temp_data,
                     file = file.path(dir_path, '00_intermidate_data', 'trace_table_clean.csv'))
  } else {
    if (length(idx_rm) == nrow(temp_data)) {
      stop(errorStyle('No features have MS2 spectrum!\n'))
    } else {
      cat(warnStyle('Warning: ', length(idx_rm), 'features without MS2 are removed.'))
    }

    dir.create(file.path(dir_path, '00_intermidate_data'), showWarnings = FALSE, recursive = TRUE)
    temp_data <- temp_data[-idx,,drop = FALSE]

    readr::write_csv(temp_data,
                     file = file.path(dir_path, '00_intermidate_data', 'trace_table_clean.csv'))
  }

}

################################################################################
# Output -----------------------------------------------------------------------
  # addEdgeAnnotation ----------------------------------------------------------
################################################################################
#' @title addEdgeLabel
#' @import stringr
addEdgeLabel <- function(tracer_table_file,
                         network_file = 'ms2_sim_table.csv',
                         id_table = id_table,
                         dir_path = '.',
                         mz_ppm = 15,
                         mass_diff_freq_cut_off = 1){

  # browser()

  tracer_table <- readr::read_csv(file.path(dir_path, '00_intermidate_data', tracer_table_file)
                                  # , show_col_types = FALSE
                                  )
  network_table <- readr::read_csv(file.path(dir_path, '00_intermidate_data', network_file)
                                   # , show_col_types = FALSE
                                   )
  id_table <- readr::read_csv(file.path(dir_path, id_table)
                              # , show_col_types = FALSE
  )

  mass_diff_table <- mass_diff_table[which(mass_diff_table$Freq >= mass_diff_freq_cut_off), ]

  # match the network_table table and id_table table to get the mz
  idx_from <- match(network_table$From_feature, id_table$name)
  network_table$from_mz <- id_table$mz[idx_from]
  idx_to <- match(network_table$To_feature, id_table$name)
  network_table$to_mz <- id_table$mz[idx_to]
  # network_table$mass_diff_abs <- abs(network_table$to_mz - network_table$from_mz)
  # network_table$mass_diff_abs_upper <- network_table$mass_diff_abs * (1 + mz_ppm * 10 ^ (-6))
  # network_table$mass_diff_abs_lower <- network_table$mass_diff_abs * (1 - mz_ppm * 10 ^ (-6))

  # the mass difference should consider the adduct forms
  idx <- str_detect(string = id_table$confidence_level, pattern = 'level1|level2') # get the level 1 and level 2 id 
  id_table_level12 <- id_table[which(idx), ]
  rm(idx)

  # idx_from <- match(network_table$From_feature, id_table_level12$peak_name)
  # idx_to <- match(network_table$To_feature, id_table_level12$peak_name)

  # the feature name of id table had been changed (20230420)
  idx_from <- match(network_table$From_feature, id_table_level12$name)
  idx_to <- match(network_table$To_feature, id_table_level12$name)
  
  network_table$temp_adduct_from <- id_table_level12$adduct[idx_from]
  network_table$temp_adduct_to <- id_table_level12$adduct[idx_to]

  if(any(stringr::str_detect(id_table_level12$adduct, '\\[M\\+H\\]\\+'))){
    network_table$temp_adduct_from[which(is.na(network_table$temp_adduct_from))] <- '[M+H]+'
    network_table$temp_adduct_to[which(is.na(network_table$temp_adduct_to))] <- '[M+H]+'
    pos <- 'positive'
  }else{
    network_table$temp_adduct_from[which(is.na(network_table$temp_adduct_from))] <- '[M-H]-'
    network_table$temp_adduct_to[which(is.na(network_table$temp_adduct_to))] <- '[M-H]-'
    pos <- 'negative'
  }


  network_table$temp_adduct_from <- sapply(seq(nrow(network_table)), function(i){
    return(strsplit(network_table$temp_adduct_from[i], ';')[[1]][1])
  })

  network_table$temp_adduct_to <- sapply(seq(nrow(network_table)), function(i){
    return(strsplit(network_table$temp_adduct_to[i], ';')[[1]][1])
  })


  if(pos == 'positive'){
    network_table$temp_adduct_from[which(network_table$temp_adduct_from == '[M]+')] <- '[M+H]+'
    network_table$temp_adduct_to[which(network_table$temp_adduct_to == '[M]+')] <- '[M+H]+'
  }else{
    network_table$temp_adduct_from[which(network_table$temp_adduct_from == '[M-2H]-')] <- '[M-H]-'
    network_table$temp_adduct_to[which(network_table$temp_adduct_to == '[M-2H]-')] <- '[M-H]-'
  }

  lib_adduct_nl <- lib_adduct_nl[[pos]]
  network_table$delta_from <- lib_adduct_nl$delta_mz[match(network_table$temp_adduct_from, lib_adduct_nl$adduct)]
  network_table$delta_to <- lib_adduct_nl$delta_mz[match(network_table$temp_adduct_to, lib_adduct_nl$adduct)]


  network_table$temp_exact_mass_from <- sapply(seq(nrow(network_table)), function(i){
    return(calculateMz(mz = network_table$from_mz[i],
                       adduct = network_table$temp_adduct_from[i],
                       delta_mz = network_table$delta_from[i]))
  }, simplify = TRUE)

  network_table$temp_exact_mass_to <- sapply(seq(nrow(network_table)), function(i){
    return(calculateMz(mz = network_table$to_mz[i],
                       adduct = network_table$temp_adduct_to[i],
                       delta_mz = network_table$delta_to[i]))
  }, simplify = TRUE)


  network_table$mass_diff_abs <- abs(as.numeric(network_table$temp_exact_mass_from) - as.numeric(network_table$temp_exact_mass_to))

  network_table$temp_adduct_from <- NULL
  network_table$temp_adduct_to <- NULL
  network_table$delta_from <- NULL
  network_table$delta_to <- NULL
  network_table$temp_exact_mass_from <- NULL
  network_table$temp_exact_mass_to <- NULL

  network_table$mass_diff_abs_upper <- as.numeric(network_table$mass_diff_abs + 0.002)
  network_table$mass_diff_abs_lower <- as.numeric(network_table$mass_diff_abs - 0.002)
  res_atom_diff <- lapply(seq(nrow(network_table)), function(i){
    temp_row <- network_table[i, ]
    temp_idx <- which(mass_diff_table$GroupDeltaMass_abs >= temp_row$mass_diff_abs_lower &
                        mass_diff_table$GroupDeltaMass_abs <= temp_row$mass_diff_abs_upper)
    # mass_accuracy <- round(abs(mass_diff_table$GroupDeltaMass_abs[temp_idx] - temp_row$mass_diff_abs) / temp_row$mass_diff_abs * 10 ^6,
    #                        digits = 0)
    DetlaMassAccuracy <- round(abs(mass_diff_table$GroupDeltaMass_abs[temp_idx] - temp_row$mass_diff_abs), digits = 4)

    # res_mass_diff_thero <- paste0(mass_diff_table$GroupDeltaMass[temp_idx],
    #                               collapse = ';')

    DetlaMassAccuracy <- paste0(DetlaMassAccuracy,
                                collapse = ';')

    AtomDiff <- paste0(mass_diff_table$AtomDifference[temp_idx],
                                 collapse = ';')

    ReactionFreq <- paste0(mass_diff_table$Freq[temp_idx],
                                 collapse = ';')

    ReactionClass <- paste0(mass_diff_table$reaction_class[temp_idx],
                                 collapse = ';')

    ReactionClassAll <- paste0(mass_diff_table$reaction_class_all[temp_idx],
                             collapse = ';')

    # with_C_diff <- paste0(mass_diff_table$with_C_diff[temp_idx],
    #                          collapse = ';')

    return(data.frame(DetlaMassAccuracy = DetlaMassAccuracy,
                      AtomDiff = AtomDiff,
                      ReactionClass = ReactionClass,
                      ReactionClassAll = ReactionClassAll,
                      ReactionFreq = ReactionFreq,
                      stringsAsFactors = FALSE))
  })

  res_atom_diff <- res_atom_diff %>% dplyr::bind_rows()
  network_table$mass_diff_abs_upper <- NULL
  network_table$mass_diff_abs_lower <- NULL
  colnames(network_table)[which(colnames(network_table) == 'mass_diff_abs')] <- 'DetlaMassAbs'

  # result <- network_table %>%
  #   dplyr::bind_cols(res_atom_diff) %>%
  #   dplyr::mutate(edge_name = paste0(From_feature," (interacts with) ", To_feature))

  result <- network_table %>%
    dplyr::bind_cols(res_atom_diff)

  # dir.create(file.path(dir_path, ), showWarnings = FALSE, recursive = TRUE)
  # readr::write_csv(result,
  #                  file = file.path(dir_path, 'edge_table_with_label.csv'))
  readr::write_csv(result,
                   file = file.path(dir_path, '00_intermidate_data', 'edge_table_with_label.csv'))

}


  # addNodeLabel ---------------------------------------------------------------

addNodeLabel <- function(tracer_table_file,
                         network_file = 'ms2_sim_table.csv',
                         dir_path = '.'){
  tracer_table <- readr::read_csv(file.path(dir_path,  '00_intermidate_data', tracer_table_file)
                                  # , show_col_types = FALSE
                                  )
  network_table <- readr::read_csv(file.path(dir_path, '00_intermidate_data', network_file)
                                   # , show_col_types = FALSE
                                   )
  node_list <- c(network_table$From_feature, network_table$To_feature) %>% unique()
  node_table <- match(node_list,
                      tracer_table$name) %>%
    tracer_table[.,] %>%
    dplyr::select(name:label)

  idx <- node_table$id %>%
    stringr::str_detect(pattern = 'NA') %>%
    which()

  node_table$attribute <- 'seed'
  node_table$attribute[idx] <- 'unknown'

  # dir.create(file.path(dir_path), showWarnings = FALSE, recursive = TRUE)
  # readr::write_csv(node_table,
  #                  file = file.path(dir_path, 'node_table_with_label.csv'))

}

################################################################################
#' @title calculateMz
#' @import stringr
calculateMz <- function(mz,
                        adduct,
                        delta_mz,
                        nmol = NULL,
                        ncharge = NULL){
  if (length(nmol) == 0) {
    if (stringr::str_detect(adduct, pattern = '2M')) {
      exact_mass <- (mz - delta_mz)/2
    } else if (stringr::str_detect(adduct, pattern = '3M')) {
      exact_mass <- (mz - delta_mz)/3
    } else {
      exact_mass <- mz - delta_mz
    }
  } else {
    exact_mass <- (mz - delta_mz)/nmol
  }

  exact_mass
}

