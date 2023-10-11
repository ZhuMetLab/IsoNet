################################################################################
# calMs2Net --------------------------------------------------------------------

#' @title calMS2Net
#' @author Zhiwei Zhou
#' @param network_file
#' @param dir_path
#' @param scoring_approach
#' @param mz_tol_ms2
#' @export

calMS2Net <- function(network_file = 'GlobalMatchResult_filterred.csv',
                      ms2_file = 'ms2_data.msp',
                      dir_path = '.',
                      scoring_approach = 'gnps',
                      mz_tol_ms2 = 15,
                      ms2_score_cutoff = 0.5,
                      num_matched_frag_cutoff = 1) {

  cat('Start calculate MS2 similarity ...\n')
  # browser()

  # read files
  pairs <- readr::read_csv(file.path(dir_path, '00_intermidate_data', network_file)
                           # , show_col_types = FALSE
                           )
  # using split by '|' instead str_replace to split feature_name and kegg_id (20230420_LMD)
  # pairs <- pairs %>%
  #   dplyr::mutate(From_feature = stringr::str_replace(From, pattern = '(NA\\d+|C\\d+|D\\d+|Kegg)', replacement = ''),
  #                 To_feature = stringr::str_replace(To, pattern = '(NA\\d+|C\\d+|D\\d+|Kegg)', replacement = ''))
  
  # pairs <- pairs %>%
  #   dplyr::mutate(From_feature = stringr::str_split(From, pattern = '\\|')[[1]][1],
  #                 To_feature = stringr::str_split(To, pattern = '\\|')[[1]][1])
  
  
  pairs$From_feature <- sapply(pairs$From, function(ff){
    stringr::str_split(ff, pattern = '\\|')[[1]][1]
  }, USE.NAMES = F)

  pairs$To_feature <- sapply(pairs$To, function(tf){
    stringr::str_split(tf, pattern = '\\|')[[1]][1]
  }, USE.NAMES = F)
  
  ms2_spec <- readMSP(file.path(dir_path, ms2_file), mode = 'all')
  names(ms2_spec) <- sapply(ms2_spec, function(x){
    x$info$NAME
  })

  # run MS/MS match
  result <- pbapply::pblapply(seq_along(pairs$Type), function(i){
    # cat(i, ' ') #2144
    name_fea1 <- pairs$From_feature[i]
    name_fea2 <- pairs$To_feature[i]
    specExp <- ms2_spec[[name_fea1]] %>% convertSpectraData()
    specRef <- ms2_spec[[name_fea2]] %>% convertSpectraData()
    score_gnps <- runSpecMatch(obj_ms2_cpd1 = specExp,
                               obj_ms2_cpd2 = specRef,
                               mz_tol_ms2 = 15,
                               scoring_approach = 'gnps')
  })

  result <- result %>% dplyr::bind_rows()
  result1 <- result %>% dplyr::select(-(name:mz))
  rownames(result1) <- NULL

  gnps_score <- pairs %>%
    dplyr::bind_cols(result1) %>%
    dplyr::rename(ms2_score = score)

  dir.create(file.path(dir_path, '00_intermidate_data'), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(gnps_score,
                   file = file.path(dir_path, '00_intermidate_data', 'ms2_sim_table.csv'))

  gnps_score <- gnps_score %>%
    dplyr::filter(ms2_score >= ms2_score_cutoff & n_frag_total >= num_matched_frag_cutoff)

  dir.create(file.path(dir_path, '00_intermidate_data'), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(gnps_score,
                   file = file.path(dir_path, '00_intermidate_data', 'ms2_sim_table_filter.csv'))

}






################################################################################
# convertSpectraData ---------------------------------------------------------
#' @title convertSpectraData
#' @param ms2_data
#' @importClassesFrom SpectraTools 'SpectraData'


convertSpectraData <- function(ms2_data) {
  options(readr.num_columns = 0)
  temp_info <- ms2_data$info %>%
    dplyr::rename(name = NAME,
                  mz = PRECURSORMZ) %>%
    dplyr::select(name:mz) %>%
    readr::type_convert()

  temp_ms2_data <- ms2_data$spec

  result <- new('SpectraData',
                info = temp_info,
                spectra = list(temp_ms2_data))

  return(result)
}


# runSpecMatch -----------------------------------------------------------------

#' @title runSpecMatch
#' @description a interphace of runing SpectraTools
#' @author Zhiwei Zhou
#' @param obj_ms2_cpd1 experimental ms2 object. Note: info is a data.frame, spec is a list formed by matrix
#' @param obj_ms2_cpd2 library ms2 object.
#' @param mz_tol_ms2 Default: 35 ppm
#' @param scoring_approach 'dp', 'bonanza', 'hybrid', 'gnps'


runSpecMatch <- function(
    obj_ms2_cpd1,
    obj_ms2_cpd2,
    mz_tol_ms2 = 35,
    scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
    ...
) {
  # browser()
  match.arg(scoring_approach)

  switch (scoring_approach,
          'dp' = {
            intensityNormedMethod <- 'maximum'
            methodScore <- 'dp'
          },
          'bonanza' = {
            intensityNormedMethod <- 'bonanza'
            methodScore <- 'bonanza'
          },
          'hybrid' = {
            intensityNormedMethod <- 'maximum'
            methodScore <- 'hybrid'
          },
          'gnps' = {
            intensityNormedMethod <- 'gnps'
            methodScore <- 'gnps'
          }
  )

  matchParam <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
                                         cutoff = 0,
                                         weightIntensity = 1,
                                         weightMZ = 0,
                                         normIntensity = TRUE,
                                         tuneLibSpectra = TRUE,
                                         intensityExpNormed = TRUE,
                                         intensityLibNormed = TRUE,
                                         includePrecursor = FALSE,
                                         ppmPrecursorFilter = 20,
                                         thrIntensityAbs = 0,
                                         thrIntensityRel = 0,
                                         intensityNormedMethod = intensityNormedMethod,
                                         methodMatch = 'direct',
                                         methodScore = methodScore) %>%
    new(Class = 'MatchParam')


  result <- try(SpectraTools::MatchSpectra(dataExp = obj_ms2_cpd1,
                                           dataRef = obj_ms2_cpd2,
                                           matchParam),
                silent = TRUE)

  if (length(result) == 0) {
    temp_result <- data.frame(name = obj_ms2_cpd2@info$name,
                              mz = obj_ms2_cpd2@info$mz,
                              score = 0,
                              n_frag_cpd1 = 0,
                              n_frag_cpd2 = 0,
                              n_frag_match = 0,
                              n_frag_nl = 0,
                              n_frag_total = 0)

    return(temp_result)
  }


  stat_matched_frag <- lapply(seq_along(result@matchedFragments), function(i){
    temp_matchedFragments <- result@matchedFragments[[i]]
    temp_nlFragments <- result@nlFragments[[i]]

    if (length(temp_matchedFragments) > 0) {
      n_frag_cpd1 <- temp_matchedFragments %>%
        dplyr::filter(intensity > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

      n_frag_cpd2 <- temp_matchedFragments %>%
        dplyr::filter(intensityExp > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

      n_frag_match <- temp_matchedFragments %>%
        dplyr::filter(intensity > 0 & intensityExp > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

    } else {
      n_frag_cpd1 <- nrow(obj_ms2_cpd1@spectra[[1]])
      n_frag_cpd2 <- nrow(obj_ms2_cpd2@spectra[[i]])
      n_frag_match <- 0
    }

    if (scoring_approach == 'dp') {
      n_nl_match <- 0
    } else {
      if (length(temp_nlFragments) > 0) {
        n_nl_match <- temp_nlFragments %>%
          dplyr::filter(intensity > 0 & intensityExp > 0) %>%
          dplyr::count() %>%
          dplyr::pull()
      } else {
        n_nl_match <- 0
      }
    }

    temp_result <- tibble::tibble(n_frag_cpd1 = n_frag_cpd1,
                                  n_frag_cpd2 = n_frag_cpd2,
                                  n_frag_match = n_frag_match,
                                  n_frag_nl = n_nl_match) %>%
      dplyr::mutate(n_frag_total = n_frag_match + n_frag_nl)

  })

  stat_matched_frag <- stat_matched_frag %>% dplyr::bind_rows()

  result@info <- result@info %>%
    dplyr::bind_cols(stat_matched_frag)

  return(result@info)
}




# readMSP ----------------------------------------------------------------------
#' @title readMSP
#' @description read MSP spectra files
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param file the file name
#' @param mode standard: extract name, mz and RT from the file, which fit for MSP data exported from QI software; all: extract all information in the file
#' @return
#' A list object. info: MS1 information, including ms1, rt etc.; spec: the MS/MS spectrum
#' @examples
#' test <- readMSP(file = 'F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/zhumetlib_validation_pos_20v_190520.msp', mode = 'all')

setGeneric('readMSP', function(file,
                               mode=c('standard', 'all')) {
  # devtools::use_package('dplyr')

  mode <- match.arg(mode)
  msp.data.list <- ListDB(file)
  nr.num.pk <- grep('Num Peaks', stringr::str_to_title(msp.data.list[[1]]))
  info.spec <- lapply(msp.data.list, function(msp.data) {
    info.list <- strsplit(msp.data[1:nr.num.pk], split = ': ')
    info <- do.call(cbind, info.list)
    colnames(info) <- info[1, ]

    if (mode=='standard') {
      mz <- round(as.numeric(info[2,3]), digits = 4)
      rt <- round(as.numeric(strsplit(info[2, "Comment"], split = "_")[[1]][1])*60, digits = 0)
      name <- info[2, "Comment"]
      # info <- matrix(c(mz, rt), ncol = 2)
      info <- data.frame(mz=mz, rt=rt)
      rownames(info) <- name
      # colnames(info) <- c("mz", "rt")
    } else {
      info <- as.data.frame(tibble::as.tibble(info))
      info <- info[-1,,drop=F]
      rownames(info) <- NULL
    }

    spec.list <- strsplit(msp.data[(nr.num.pk+1):length(msp.data)], ' |\\t')
    spec.list <- lapply(spec.list, as.numeric)
    spec <- do.call(rbind, spec.list)[, c(1, 2), drop = FALSE]
    colnames(spec) <- c('mz', 'intensity')
    # spec <- list('spec' = spec)

    # list('info' = info[-1, , drop = FALSE],
    #      'spec' = spec)

    list('info' = info,
         'spec' = spec)

  })

  return(info.spec)
})

setGeneric('ListDB', function(file) {
  msp.data <- readLines(file)
  nl.db.new <- 1
  idx.db <- 1
  db.list <- list()
  len.data <- length(msp.data)
  for(nl in 1:len.data)
  {
    if(msp.data[nl]=="") {
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:(nl-1)])
      nl.db.new <- nl + 1
      idx.db <- idx.db + 1
    } else if (nl == len.data){
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:nl])
    }
  }
  db.list
})

