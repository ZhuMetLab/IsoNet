################################################################################
# calMIDSim --------------------------------------------------------------------

#' @title calMIDSim
#' @author Hongmiao Wang, Zhiwei Zhou
#' @param NodeTable
#' @param Pattern_name
#' @param TargetTable_dir Stable isotope labeling quantification table. CSV table. First 1-6 column: name, id, mz, compound, label. The label column: 0-n represent isotopologue M0 - Mn
#' @param Cutoff MID similarity score. Default: 0.9
#' @param IsoDegree minimum labeling extent. Note: the labeled isotopologue can't be reduced in similarity scoring. Default: 0.1
#' @param FoldChange a parameter to judge similarity scoring type. If the most abundant isotopologue divide the second abundant isotopologue larger than this value, the isotopologue simialrity is type I, else the type is type II. Default: 20
#' @param MinMotifLength Character. Minimum length of motif (value < 1 means fraction of max motif length). Default: 0.5
#' @param maxMotif the maximum motif combination. Default: 10000
#' @return
#' @export

# mid_cutoff = 0.95
# mid_fc = 20
# mid_isoDegree = 0.02
# mid_min_motifLen = 0.5
# min_max_motif = 1000

# NodeTable = NULL
# Pattern_name = "Unknown"
# TargetTable_dir = file.path('00_intermidate_data', 'trace_table_clean.csv')
# TargetSample = c("GGA_1","GGA_2","GGA_3","GGA_4","GGA_5","GGA_6")
# dir_path = 'G:/00_projects/06_isotopeTracing/data_processing/20220831_MIDNet'
# Cutoff = 0.9
# IsoDegree = 0.1
# FoldChange = 20
# MinMotifLength = 0.5
# multi_thread = 1
# maxMotif = 10000


calMIDSim <- function(NodeTable = NULL,
                      Pattern_name = "Unknown",
                      TargetTable_dir = "./QuantTableContiminateAdjustScaled.csv",
                      TargetSample = NULL,
                      dir_path = '.',
                      Cutoff = 0.9,
                      IsoDegree = 0.1,
                      FoldChange = 20,
                      MinMotifLength = 0.5,
                      multi_thread = 1,
                      maxMotif = 10000,
                      ignore_max = FALSE)
{
  # browser()
  # generate possible pairs
  TargetTable <- read.csv(file.path(dir_path, TargetTable_dir))
  TargetTable$id[which(is.na(TargetTable$id))] <- 'NA1'
  TargetTable[is.na(TargetTable)] <- 0
  TargetTable$name_match <- paste0(TargetTable$name, '|', TargetTable$id ) # modified (previous: no '|' split sign)
  featurelist <- as.character(unique(TargetTable$name_match))
  # featurelist <- unique(as.character(TargetTable$name_match))
  if(is.null(NodeTable))
  {
    featureComb <- combn(featurelist,2)
    NodeTable <- as.data.frame(t(featureComb))
    colnames(NodeTable) <- c("From","To")
  }
  TypeList <- c()
  C_number <- c()
  for (i in seq_along(featurelist)) {
    FeatureTable <- TargetTable[which(TargetTable$name_match==featurelist[i]),]
    QuantTable <- FeatureTable[,TargetSample]
    QuantTable <- QuantTable[-1,]
    QuantValue <- as.numeric(apply(QuantTable, 1, mean))
    Rank <- order(QuantValue, decreasing=TRUE)
    C_number <- c(C_number,length(QuantValue))
    if(QuantValue[Rank[1]]!=0)
    {
      if(QuantValue[Rank[1]]/QuantValue[Rank[2]] > FoldChange)
      {
        TypeList <- c(TypeList,"TypeII")
      }else{
        TypeList <- c(TypeList,"TypeIII")
      }
    }else{
      TypeList <- c(TypeList,"non-Type")
    }
  }
  # NodeTable$Type <- "non-Type"

  # add types and only keep same type pairs
  cat("Check NodeTable. \n")
  idx_from = match(NodeTable$From, featurelist)
  idx_to = match(NodeTable$To, featurelist)
  NodeTable <- NodeTable %>%
    dplyr::mutate(from_type = TypeList[idx_from],
                  to_type = TypeList[idx_to],from_C = C_number[idx_from],to_C = C_number[idx_to]) %>%
    dplyr::mutate(type_label = (from_type == to_type)|(from_C==to_C)) %>%
    dplyr::filter(type_label)%>%
    dplyr::mutate(Type = from_type)
  NodeTable$C_diff <- NodeTable$from_C-NodeTable$to_C
  NodeTable$Type[which(NodeTable$C_diff==0)] <- "TypeI"
  NodeTable <- NodeTable %>% dplyr::select(-c('from_type', 'to_type', 'type_label', 'from_C', 'to_C','C_diff'))


  # MaxScoreAll <- c()
  # BestMotif <- c()
  # BestMatch <- c()
  # MotifLength <- c()
  # MotifSource <- c()
  cat("Matching Start. \n")
  start_time <- Sys.time()
  cat("Calculate Scores. \n")
  # pb <- txtProgressBar(style = 3)

  # BiocParallel::bpmapply(function(), )

  if (multi_thread == 1) {
    list_SIMres <- pbapply::pbmapply(function(patternA, patternB, type){
      SIMResAll <- MIDSIM(patternA,
                          patternB,
                          TargetTable = TargetTable,
                          TargetSample = TargetSample,
                          Cutoff = Cutoff,
                          IsoDegree = IsoDegree,
                          MinMotifLength = MinMotifLength,
                          printall = FALSE,
                          printCos = FALSE,
                          type = type,
                          ignore_max=ignore_max)
    },
    patternA = as.character(NodeTable$From),
    patternB = as.character(NodeTable$To),
    type = as.character(NodeTable$Type),
    SIMPLIFY = FALSE)

    # list_SIMres1 <- pbapply::pbmapply(function(patternA, patternB, type){
    #   SIMResAll <- MIDSIM(patternA,
    #                       patternB,
    #                       TargetTable = TargetTable,
    #                       TargetSample = TargetSample,
    #                       Cutoff = Cutoff,
    #                       IsoDegree = IsoDegree,
    #                       MinMotifLength = MinMotifLength,
    #                       printall = FALSE,
    #                       printCos = FALSE,
    #                       type = type)
    # },
    # patternA = as.character(NodeTable$From)[1:1000],
    # patternB = as.character(NodeTable$To)[1:1000],
    # type = as.character(NodeTable$Type)[1:1000],
    # SIMPLIFY = FALSE)

  } else {

    # Parallel run mapply
    # library(parallel)
    # cl <- makeCluster(multi_thread)
    #
    # tempFun <- function(patternA, patternB, type) {
    #   SIMResAll <- MIDSIM(patternA,
    #                       patternB,
    #                       TargetTable = TargetTable,
    #                       TargetSample = TargetSample,
    #                       Cutoff = Cutoff,
    #                       IsoDegree = IsoDegree,
    #                       MinMotifLength = MinMotifLength,
    #                       printall = FALSE,
    #                       printCos = FALSE,
    #                       type = type)
    # }
    #
    # clusterExport(cl,
    #               c("tempFun", "MIDSIM", 'NodeTable', 'TargetTable', 'TargetSample', 'Cutoff', 'IsoDegree', 'MinMotifLength', '%>%'),
    #               envir = environment())
    #
    # system.time(
    #   list_SIMres <- mcmapply(tempFun,
    #                           patternA = as.character(NodeTable$From),
    #                           patternB = as.character(NodeTable$To),
    #                           type = as.character(NodeTable$Type),
    #                           # mc.cores = 4L,
    #                           SIMPLIFY = FALSE)
    # )
    #
    # stopCluster(cl)


    # BiocParallel run mapply
    tempFun <- function(i, patternA, patternB, type, ignore_max, TargetTable, TargetSample, Cutoff, IsoDegree, MinMotifLength,
                        MIDSIM, DataInput, MIDSimilarityCal, ConvertMotif, CalCosine, calMotifLength) {
      temp_patternA <- patternA[i]
      temp_patternB <- patternB[i]
      temp_type <- type[i]

      SIMResAll <- MIDSIM(temp_patternA,
                          temp_patternB,
                          TargetTable = TargetTable,
                          TargetSample = TargetSample,
                          Cutoff = Cutoff,
                          IsoDegree = IsoDegree,
                          MinMotifLength = MinMotifLength,
                          printall = FALSE,
                          printCos = FALSE,
                          type = temp_type,
                          ignore_max = ignore_max)
    }

    snow <- BiocParallel::SnowParam(workers = multi_thread,
                                    type = "SOCK",
                                    progressbar = TRUE)

    list_SIMres <- BiocParallel::bplapply(seq_along(NodeTable$From), # 1:10000
                                          FUN = tempFun,
                                          BPPARAM = snow,
                                          patternA = as.character(NodeTable$From),
                                          patternB = as.character(NodeTable$To),
                                          type = as.character(NodeTable$Type),
                                          ignore_max=ignore_max,
                                          TargetTable = TargetTable,
                                          TargetSample = TargetSample,
                                          Cutoff = Cutoff,
                                          IsoDegree = IsoDegree,
                                          MinMotifLength = MinMotifLength,
                                          MIDSIM = MIDSIM,
                                          DataInput = DataInput,
                                          MIDSimilarityCal = MIDSimilarityCal,
                                          ConvertMotif = ConvertMotif,
                                          CalCosine = CalCosine,
                                          calMotifLength = calMotifLength)
  }




#
#
#
#   end_time <- Sys.time()
#   cat("Matching End. Cost ", end_time - start_time, "\n")
#
  # list_SIMres <- pbapply::pblapply(seq(10000), function(i){
  #   patternA <- as.character(NodeTable[i,1])
  #   patternB <- as.character(NodeTable[i,2])
  #   SIMResAll <- MIDSIM(patternA,
  #                       patternB,
  #                       TargetTable = TargetTable,
  #                       TargetSample = TargetSample,
  #                       Cutoff = Cutoff,
  #                       IsoDegree = IsoDegree,
  #                       MinMotifLength = MinMotifLength,
  #                       printall = FALSE,
  #                       printCos = FALSE,
  #                       type = as.character(NodeTable[i,3]))
  # })
  #




  # retrieve motif source
  SIMRes <- lapply(list_SIMres, function(x){x[[1]]})
  POrder <- sapply(list_SIMres, function(x){x[[2]]})
  type <- sapply(list_SIMres, function(x){x[[3]]})

  temp_idx <- which(POrder == 2)
  MotifSource <- NodeTable$To
  MotifSource[temp_idx] <- NodeTable$From[temp_idx]
  Reference <- NodeTable$From
  Reference[temp_idx] <- NodeTable$To[temp_idx]
  
  MaxScoreAll <- sapply(SIMRes, function(x){x[[1]]})
  BestMotif <- sapply(SIMRes, function(x){x[[3]]})

  MotifLength <- sapply(BestMotif, function(x){
    if (grepl("-",x)) {
      MotifLength <- calMotifLength(x)
    } else {
      MotifLength <- 0
    }
    return(MotifLength)
  })
  
  BestMatch <- sapply(SIMRes, function(x){x[[2]]})
  ReferenceRange <- paste0('M',BestMatch,"-M", (BestMatch+MotifLength-1))
  
  #BestMatch <- paste0('M', sapply(SIMRes, function(x){x[[2]]}))
  
  

  NodeTable$Type <- type
  NodeTable <- NodeTable %>%
    dplyr::mutate(MIDSimilarityScore = MaxScoreAll,
                  MotifSource = MotifSource,
                  BestMotif = BestMotif,
                  MotifLength = MotifLength,
                  Reference = Reference,
                  ReferenceRange = ReferenceRange)

  dir.create(file.path(dir_path, '00_intermidate_data'), showWarnings = FALSE, recursive = TRUE)
  save(list_SIMres, file = file.path(dir_path, '00_intermidate_data', 'list_SIMres.RData'))

  readr::write_csv(NodeTable,
                   file = file.path(dir_path, '00_intermidate_data', "GlobalMatchResult.csv"))

  NodeTable <- NodeTable %>%
    dplyr::filter(MIDSimilarityScore >= Cutoff)
  readr::write_csv(NodeTable,
                   file = file.path(dir_path, '00_intermidate_data', "GlobalMatchResult_filterred.csv"))

  # # for (i in 1:nrow(NodeTable)) {
  # for (i in 1:10000) {
  #   # i 13265
  #   # cat(i, ' ')
  #   # setTxtProgressBar(pb, i/nrow(NodeTable))
  #   setTxtProgressBar(pb, i/10000)
  #   patternA <- as.character(NodeTable[i,1])
  #   patternB <- as.character(NodeTable[i,2])
  #   SIMResAll <- MIDSIM(patternA,
  #                       patternB,
  #                       TargetTable = TargetTable,
  #                       TargetSample = TargetSample,
  #                       Cutoff = Cutoff,
  #                       IsoDegree = IsoDegree,
  #                       MinMotifLength = MinMotifLength,
  #                       printall = FALSE,
  #                       printCos = FALSE,
  #                       type = as.character(NodeTable[i,3]))
  #   SIMRes <- SIMResAll[[1]]
  #   POrder <- SIMResAll[[2]]
  #   if (POrder==1) {
  #     MotifSource <- c(MotifSource, as.character(NodeTable[i,2]))
  #   } else {
  #     MotifSource <- c(MotifSource, as.character(NodeTable[i,1]))
  #   }
  #   MaxScoreAll <- c(MaxScoreAll, SIMRes[[1]])
  #   if (length(SIMRes)==0 | length(SIMRes[[3]])==0) {
  #     BestMotif <- c(BestMotif,"No Matched")
  #   } else {
  #     BestMotif <- c(BestMotif,SIMRes[[3]])
  #   }
  #   if (length(SIMRes[[3]])==0)
  #     SIMRes[[3]] <- "NA"
  #   if (grepl("-",SIMRes[[3]])) {
  #     MotifLength <- c(MotifLength,calMotifLength(SIMRes[[3]]))
  #   } else {
  #     MotifLength <- c(MotifLength,0)
  #   }
  #   BestMatch <- c(BestMatch, paste0("M",SIMRes[[2]]))
  # }
  #
  # end_time <- Sys.time()
  # close(pb)
  # cat("Matching End. Cost ",end_time-start_time,"\n")
  # NodeTable$MaxScore <- MaxScoreAll
  # NodeTable$BestMotif <- BestMotif
  # NodeTable$MotifLength <- MotifLength
  # NodeTable$BestMatch <- BestMatch
  # NodeTable$MotifSource <- MotifSource
  # write.csv(NodeTable,"GlobalMatchResult.csv")
  # NodeTable <- dplyr::filter(NodeTable, MaxScore>=Cutoff)
  # write.csv(NodeTable,"GlobalMatchResult_filterred.csv")

}

################################################################################
  # MIDSIM ---------------------------------------------------------------------


MIDSIM <- function(patternA,
                   patternB,
                   Pattern_name = "Unknown",
                   UseShort = TRUE,
                   TargetTable = NULL,
                   TargetSample = NULL,
                   Cutoff = 0.9,
                   IsoDegree = 0.1,
                   MinMotifLength = 3,
                   printall = FALSE,
                   printCos = TRUE,
                   range = NULL,
                   GenMotif = TRUE,
                   type = "TypeIII",
                   maxMotif = 10000,
                   ignore_max = FALSE)
{
  # browser()
  library(dplyr)
  library(stringr)
  # library(MIDNet)

  if(!is.null(TargetTable))
  {
    res <- DataInput(patternA, patternB, TargetTable, TargetSample)
    patternA <- res[[1]]
    patternB <- res[[2]]
    Pattern_name <- res[[3]]
  }
  if(UseShort)
  {
    # use short one as motif
    if(length(patternA)>length(patternB))
    {
      PMatch <- patternB
      PSearch <- patternA
      POrder <- 1
    } else {
      PMatch <- patternA
      PSearch <- patternB
      POrder <- 2
    }
  } else {
    PMatch <- patternA
    PSearch <- patternB
    POrder <- 2
  }
  # if(length(patternA)==length(patternB))
  # {
  #   type="TypeI"
  # }

  if(type=="TypeI"|type=="TypeII")
  {
    if(which.max(PSearch[-1]) %in% 1:length(PMatch[-1]))
    {

      if(ignore_max&(length(PSearch)==length(PMatch)))
      {
        GenMotif <- FALSE
      }else{

        if(which.max(PSearch[-1])==which.max(PMatch[-1]))
        {
          PSearch <- PSearch[1:length(PMatch)]
          GenMotif <- FALSE
        }else{
          SIMRes <- list(0,0,"Unmatched type I Motif in carbon index")
          return(list(SIMRes, POrder, type))
        }

      }
    }else{
      SIMRes <- list(0,0,"Unmatched type I Motif")
      return(list(SIMRes, POrder, type))
    }
  }else{
    # remove M0
    PSearch <- PSearch[-1]
    PMatch <- PMatch[-1]
  }
  if(!is.null(range))
  {
    PMatch <- PMatch[range[1]:range[2]]
    if(length(PMatch)>length(PSearch))
    {SIMRes <- list(0,0,"Longer Match Motif")
    return(list(SIMRes,POrder, type))}
  }
  SIMRes <- MIDSimilarityCal(Pattern_name,
                             PSearch,
                             PMatch,
                             Cutoff,
                             IsoDegree,
                             MinMotifLength,
                             printall,
                             GenMotif = GenMotif,
                             maxMotif = maxMotif)
  if (printCos)
  {
    cat("The match result: \n")
    cat("-------------- \n")
    cat("MaxScore:", SIMRes[[1]],"\n")
    cat("Best Motif:", SIMRes[[3]],"\n")
    cat("Best Match: ", paste0("M",SIMRes[[2]]), "\n")
    cat("-------------- \n")
  } else {
    return(list(SIMRes, POrder, type))
  }
}

  # DataInput ------------------------------------------------------------------

DataInput <- function(patternA,
                      patternB,
                      TargetTable = NULL,
                      TargetSample = NULL)
{
  if(is.null(TargetTable$name_match))
    TargetTable$name_match <- paste0(TargetTable$name,TargetTable$id)
  Pattern_name <- paste0(patternA,"_",patternB)
  TargetTable$name_match <- as.character(TargetTable$name_match)
  patternATable <- dplyr::filter(TargetTable,name_match==patternA)
  patternATable <- patternATable[,TargetSample]
  patternA <- apply(patternATable, 1, mean)
  patternBTable <- dplyr::filter(TargetTable,name_match==patternB)
  patternBTable <- patternBTable[,TargetSample]
  patternB <- apply(patternBTable, 1, mean)
  return(list(patternA,patternB,Pattern_name))
}




  # MIDSimilarityCal -----------------------------------------------------------
MIDSimilarityCal <- function(Pattern_name,
                             PSearch,
                             PMatch,
                             Cutoff = 0.9,
                             IsoDegree = 0.1,
                             MinMotifLength = 3,
                             printall = FALSE,
                             GenMotif = TRUE,
                             maxMotif = 10000)
{
  # browser()
  if(MinMotifLength < 1)
  {
    MinMotifLength <- ceiling(MinMotifLength*max(length(PSearch),length(PMatch)))
  }
  if(length(PMatch)<MinMotifLength)
    return(list(0,0,"Low Carbon number"))
  AllMotif_res <- ConvertMotif(PMatch,
                               GenMotif=GenMotif,
                               maxMotif=maxMotif,
                               MinMotifLength=MinMotifLength,
                               IsoDegree=IsoDegree)
  MotifTable <- AllMotif_res[[1]]
  MotifString <- AllMotif_res[[2]]

  if(max(MotifTable$MotifLength)<MinMotifLength)
    return(list(0,0,"Low Carbon number"))

  # match motif with one step isotopelogue
    # 1. generate possible pairs (1 step movement)
    # 2. calculate cosine score
    # 3. export largest score

  MotifTable <- MotifTable %>% dplyr::filter(MotifLength>=MinMotifLength)

  # generate possible pairs
  len_match <- sapply(MotifString, function(x){length(PSearch)-length(x)})
  len_match[len_match == 0] <- 1
  motif_index <- MotifTable$MotifIndex
  idx_search <- mapply(function(x, y){
    rep(x, each=y)
  },
  x = motif_index,
  y = len_match,
  SIMPLIFY = TRUE) %>%
    unlist()
  idx_match <- sapply(len_match, function(x){seq(x)}) %>% unlist()

  # calculate score for all pairs
  score <- mapply(function(x, y){
    Motif <- MotifString[[x]]
    if (!all(Motif<IsoDegree)) {
      Score <- round(CalCosine(Motif,PSearch[y:(y+length(Motif)-1)]),4)
    } else {
      return(0)
    }
  },
  x = idx_search,
  y = idx_match,
  SIMPLIFY = TRUE)

  score[is.nan(score)] <- 0
  idx_max <- which.max(score)
  MaxScore <- score[idx_max]
  if (MaxScore > 0) {
    BestMatch <- idx_match[idx_max]
    BestMotif <- MotifTable$MotifLabel[idx_search[idx_max]]
  } else {
    BestMatch <- 0
    BestMotif <- "Low Label Degree"
  }

  result <- list(MaxScore, BestMatch, BestMotif)
  return(result)


  # MotifTable$Score = 0
  # MotifTable$MatchPlace = 0
  # MaxScore <- 0
  # BestMatch <- 0
  # BestMotif <- 0
  #
  # for(MotifID in 1:nrow(MotifTable))
  # {
  #   Motif <- MotifString[[MotifID]]
  #   if(!all(Motif<IsoDegree))
  #   {
  #     MatchTime <- length(PSearch)-length(Motif)
  #     if(MatchTime==0)
  #       MatchTime=1
  #     for (i in 1:MatchTime) {
  #       Score <- round(CalCosine(Motif,PSearch[i:(i+length(Motif)-1)]),4)
  #       if(is.na(Score))
  #         Score <- 0
  #       if(i==1)
  #       {
  #         ScoreAll <- Score
  #         MatchPlace <- i
  #       }else{
  #         ScoreAll <- paste0(ScoreAll,";",Score)
  #         MatchPlace <- paste0(MatchPlace,";",i)
  #       }
  #       if(MaxScore<Score)
  #       {
  #         MaxScore <- Score
  #         BestMatch <- i
  #         BestMotif <- MotifID
  #       }
  #     }
  #     if((MotifID==nrow(MotifTable))|(MotifTable$MotifLength[MotifID+1]!=MotifTable$MotifLength[MotifID]))
  #     {
  #       if(MaxScore>Cutoff)
  #         if(!(printall))
  #           return(list(MaxScore,BestMatch,as.character(MotifTable$MotifLabel[BestMotif])))
  #     }
  #     MotifTable$Score[MotifID] <- ScoreAll
  #     MotifTable$MatchPlace[MotifID] <- MatchPlace
  #   } else{
  #     if(!dir.exists("SIMCheck"))
  #       dir.create("SIMCheck")
  #     write.csv(MotifTable,paste0("SIMCheck/",Pattern_name,".csv"))
  #     return(list(0,0,"Low Label Degree"))
  #   }
  # }
  # if(!dir.exists("SIMCheck"))
  #   dir.create("SIMCheck")
  # write.csv(MotifTable,paste0("SIMCheck/",Pattern_name,".csv"))
  # return(list(MaxScore,BestMatch,as.character(MotifTable$MotifLabel[BestMotif])))
}



  # ConvertMotif ---------------------------------------------------------------
ConvertMotif <- function(PMatch,
                         GenMotif = TRUE,
                         maxMotif = 10000,
                         MinMotifLength = 0,
                         IsoDegree = 0.1)
{
  MotifCos <- c()
  MotifLength <- c()
  MotifLabel <- c()
  # AllMotif <- lapply(1:length(PMatch),function(i) combn(PMatch,i))
  MotifString <- list()
  LabelPool <- list(paste0("M", 1:length(PMatch)))
  StringPool <- list(PMatch)
  if(GenMotif)
  {
    LabelContain <- LabelPool[[1]][which(StringPool[[1]] > IsoDegree)]
    for(i in 1:(length(PMatch)-1))
    {
      NewStringPool <- c()
      NewLabelPool <- c()
      for(k in 1:length(StringPool))
      {
        Motif <- StringPool[[k]]
        if(length(Motif)<MinMotifLength)
          break
        Label <- LabelPool[[k]]
        MotifString <- append(MotifString,list(Motif))
        MotifCos <- c(MotifCos,paste(list(Motif),sep = ","))
        MotifLabel <- c(MotifLabel,paste0(Label[1],"-",Label[length(Label)]))
        MotifLength <- c(MotifLength,length(Motif))
        NewStringPool <- append(NewStringPool,list(Motif[-1]))
        NewStringPool <- append(NewStringPool,list(Motif[-length(Motif)]))
        NewLabelPool <- append(NewLabelPool,list(Label[-1]))
        NewLabelPool <- append(NewLabelPool,list(Label[-length(Label)]))
        if(length(MotifLength)>=maxMotif)
          break
        if(length(Motif)<MinMotifLength)
          break
      }
      if(is.null(NewLabelPool))
        break
      StringPool <- c()
      LabelPool <- c()
      for(j in 1:length(NewStringPool))
      {
        if(all(LabelContain %in% NewLabelPool[[j]]))
        {
          StringPool <- append(StringPool,list(NewStringPool[[j]]))
          LabelPool <- append(LabelPool,list(NewLabelPool[[j]]))
        }
      }
      if(is.null(LabelPool))
        break
      if(length(MotifLength)>=maxMotif)
        break
      if(length(Motif)<MinMotifLength)
        break
    }
  }else{
    Motif <- PMatch
    Label <- paste0("M", 0:(length(PMatch)-1))
    MotifString <- append(MotifString,list(Motif))
    MotifCos <- c(MotifCos,paste(list(Motif),sep = ","))
    MotifLabel <- c(MotifLabel,paste0(Label[1],"-",Label[length(Label)]))
    MotifLength <- c(MotifLength,length(Motif))
  }

  MotifIndex <- 1:length(MotifCos)
  MotifTable <- cbind.data.frame(MotifIndex,MotifCos,MotifLength,MotifLabel)
  return(list(MotifTable,MotifString))
}



  # calMotifLength -------------------------------------------------------------
calMotifLength <- function(MotifLabel)
{
  MotifLabel <- gsub("M","",MotifLabel)
  MotifLabel <- strsplit(MotifLabel,"-")
  MotifLength <- as.numeric(MotifLabel[[1]][2])-as.numeric(MotifLabel[[1]][1])+1
  return(MotifLength)
}
#   # CalCosine ------------------------------------------------------------------
# CalCosine <- function(Vector1,Vector2)
# {
#   Score <- sum(Vector1 * Vector2) / (sqrt(sum(Vector1 ^ 2)) * sqrt(sum(Vector2 ^ 2)))
#   return(Score)
# }

# CalManhattan ------------------------------------------------------------------
CalCosine <- function(Vector1,Vector2)
{
  Score <- 1/(1+sum(abs(Vector1-Vector2)))
  return(Score)
}

