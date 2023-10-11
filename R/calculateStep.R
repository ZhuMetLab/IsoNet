################################################################################
#' @title calculateStep
#' @author Mingdu Luo
#' @import igraph
#' @import stringr


calculateStep <- function(id_table, 
                          dir_path = '.'){

  
  # get the pool of feature-metabolite pairs #####################################
  id_table <- readr::read_csv(file.path(dir_path, id_table)
                              # , show_col_types = FALSE
  )
  # id_table <- readr::read_csv('E:/07_MIDNet/20230201_293T_pos/20230130_update_node_table_pos_integrated_RP_add_remove_input.csv')
  idx <- str_detect(string = id_table$confidence_level, pattern = 'level1|level2') # get the level 1 and level 2 id 
  id_table_level12 <- id_table[which(idx), ]
  # id_table_level12_filter_na <- id_table_level12[!is.na(id_table_level12$id_kegg), ]
  id_table_level12_filter_na <- id_table_level12
  id_pool <- lapply(seq(nrow(id_table_level12_filter_na)), function(i){
    # i <- 1
    temp_row <- id_table_level12_filter_na[i, ]
    # temp_feature <- temp_row$peak_name
    # temp_id <- temp_row$id_kegg
    temp_feature <- temp_row$name
    temp_id <- temp_row$id
    
    # if(str_detect(temp_id, pattern = '/')){
    #   all_id <- str_split(temp_id, '/')[[1]]
    # }else{
    #   all_id <- str_split(temp_id, ';')[[1]]
    # }
    all_id <- str_split(temp_id, ';')[[1]]
    all_formula <- str_split(temp_row$formula, ';')[[1]]
    # all_name <- str_split(temp_row$name, ';')[[1]]
    all_name <- str_split(temp_row$compound.name, ';')[[1]]
    # all_adduct <- str_split(temp_row$adduct, ';')[[1]]
    
    
    df <- data.frame(name = temp_feature, 
                     all_id = all_id,
                     all_name = all_name, 
                     all_formula = all_formula)
    return(df)
  })
  id_pool <- do.call(rbind, id_pool)
  id_pool$seq <- paste0('#', seq(nrow(id_pool)))
  # remove NA
  # id_pool <- id_pool[!id_pool$all_id == 'NA', ]
  all_feature_in_pool <- unique(id_pool$name)
  # id_pool is the final pool
  
  # assign id to the edge table ##################################################
  edge_table <- readr::read_csv(file.path(dir_path, '00_intermidate_data', 'edge_table_with_label.csv')
                                # , show_col_types = FALSE
  )
  # edge_table <- readr::read_csv('E:/07_MIDNet/20230201_293T_pos/edge_table_with_label.csv')
  
  idx_remain_edge <- edge_table$From_feature %in% all_feature_in_pool & edge_table$To_feature %in% all_feature_in_pool
  edge_table_filter <- edge_table[idx_remain_edge, ]
  
  step_result <- lapply(seq(nrow(edge_table_filter)), function(t){
    # t <- 1
    temp_from <- edge_table_filter$From_feature[t]
    temp_to <- edge_table_filter$To_feature[t]
    
    # from_ids <- id_pool$all_id[id_pool$peak_name == temp_from]
    # to_ids <- id_pool$all_id[id_pool$peak_name == temp_to]
    # 
    # all_combn <- lapply(seq_along(from_ids), function(k){
    #   return(data.frame(from_id = from_ids[k], 
    #                     to_id = to_ids))
    # })
    # all_combn <- do.call(rbind, all_combn)
    # all_combn$From_feature <- temp_from
    # all_combn$To_feature <- temp_to
    # all_combn$step <- sapply(seq(nrow(all_combn)), function(zz){
    #   # zz <- 1
    #   return(unname(distances(mrn, v = V(mrn)[name == all_combn$from_id[zz]], 
    #                           to = V(mrn)[name == all_combn$to_id[zz]])))
    # }, simplify = T)
    # 
    # all_combn$from_name <- cpd_emrn$name[match(all_combn$from_id,cpd_emrn$id)]
    # all_combn$to_name <- cpd_emrn$name[match(all_combn$to_id,cpd_emrn$id)]
    # 
    # all_combn$from_formula <- cpd_emrn$formula[match(all_combn$from_id,cpd_emrn$id)]
    # all_combn$to_formula <- cpd_emrn$formula[match(all_combn$to_id,cpd_emrn$id)]
    # all_combn$reaction_type <- 'unknown'
    # all_combn$step <- as.numeric(all_combn$step)
    # all_combn$reaction_type[which(!is.na(all_combn$step))] <- 'known'
    
    from_seq <- id_pool$seq[id_pool$name == temp_from]
    to_seq <- id_pool$seq[id_pool$name == temp_to]
    
    all_combn <- lapply(seq_along(from_seq), function(k){
      return(data.frame(from_seq = from_seq[k],
                        to_seq = to_seq))
    })
    
    all_combn <- do.call(rbind, all_combn)
    all_combn$From_feature <- temp_from
    all_combn$To_feature <- temp_to
    
    all_combn$from_name <- id_pool$all_name[match(from_seq, id_pool$seq)]
    all_combn$to_name <- id_pool$all_name[match(to_seq, id_pool$seq)]
    
    all_combn$from_id <- id_pool$all_id[match(from_seq, id_pool$seq)]
    all_combn$to_id <- id_pool$all_id[match(to_seq, id_pool$seq)]
    
    all_combn$from_formula <- id_pool$all_formula[match(from_seq, id_pool$seq)]
    all_combn$to_formula <- id_pool$all_formula[match(to_seq, id_pool$seq)]
    
    all_combn$step <- sapply(seq(nrow(all_combn)), function(zz){
      # browser()
      res <- try(unname(distances(mrn, v = V(mrn)[name == all_combn$from_id[zz]],
                                  to = V(mrn)[name == all_combn$to_id[zz]])), 
                silent = TRUE)
      if(class(res)[1] == 'try-error'){
        return(NA)
      }else{
        return(res)
      }
    }, simplify = T)
    
    # all_combn$reaction_type <- 'unknown_1'
    all_combn$reaction_type <- 'Known'
    # all_combn$step <- as.numeric(all_combn$step)
    # all_combn$reaction_type[which(!is.na(all_combn$step))] <- 'known'
    
    return(all_combn)
  })
  
  step_result_combined <- lapply(step_result, function(res){
    # browser()
    res_combiend <- data.frame(From_feature = res$From_feature[1], 
                               To_feature = res$To_feature[1],
                               from_id = paste0(res$from_id, collapse = ';'), 
                               to_id = paste0(res$to_id, collapse = ';'), 
                               step = paste0(res$step, collapse = ';'), 
                               from_name = paste0(res$from_name, collapse = ';'), 
                               to_name = paste0(res$to_name, collapse = ';'), 
                               from_formula = paste0(res$from_formula, collapse = ';'), 
                               to_formula = paste0(res$to_formula, collapse = ';'), 
                               reaction_type = paste0(res$reaction_type, collapse = ';'))
    return(res_combiend)
  })
  
  step_result_combined <- do.call(rbind, step_result_combined)
  step_result_combined$from_to_feature <- paste0(step_result_combined$From_feature, 
                                                 '_', 
                                                 step_result_combined$To_feature)
  edge_table$from_to_feature <- paste0(edge_table$From_feature, 
                                       '_', 
                                       edge_table$To_feature)
  
  edge_table_add_id <- dplyr::left_join(edge_table, 
                                        step_result_combined, 
                                        by = 'from_to_feature')
  edge_table_add_id$from_to_feature <- NULL
  edge_table_add_id$From_feature.y <- NULL
  edge_table_add_id$To_feature.y <- NULL
  
  # edge_table_add_id$reaction_type[which(is.na(edge_table_add_id$reaction_type))] <- 'unknown_3'
  edge_table_add_id$reaction_type[which(is.na(edge_table_add_id$reaction_type))] <- 'Unknown2'
  colnames(edge_table_add_id)[which(colnames(edge_table_add_id) == 'From_feature.x')] <- 'From_feature'
  colnames(edge_table_add_id)[which(colnames(edge_table_add_id) == 'To_feature.x')] <- 'To_feature'
  idx_1 <- which(edge_table_add_id$From_feature %in% all_feature_in_pool & !edge_table_add_id$To_feature %in% all_feature_in_pool)
  idx_2 <- which(!edge_table_add_id$From_feature %in% all_feature_in_pool & edge_table_add_id$To_feature %in% all_feature_in_pool)
  
  idx_1_1 <- match(edge_table_add_id$From_feature[idx_1], id_pool$name)
  edge_table_add_id$from_id[idx_1] <- id_pool$all_id[idx_1_1]
  edge_table_add_id$from_name[idx_1] <- id_pool$all_name[idx_1_1]
  edge_table_add_id$from_formula[idx_1] <- id_pool$all_formula[idx_1_1]
  
  idx_2_2 <- match(edge_table_add_id$To_feature[idx_2], id_pool$name)
  edge_table_add_id$to_id[idx_2] <- id_pool$all_id[idx_2_2]
  edge_table_add_id$to_name[idx_2] <- id_pool$all_name[idx_2_2]
  edge_table_add_id$to_formula[idx_2] <- id_pool$all_formula[idx_2_2]
  
  
  all_idx <- c(idx_1, 
               idx_2)
  
  # edge_table_add_id$reaction_type[all_idx] <- 'unknown_2'
  edge_table_add_id$reaction_type[all_idx] <- 'Unknown1'
  
  readr::write_csv(edge_table_add_id,
                   file = file.path(dir_path, '00_intermidate_data', 'edge_table_with_label_add_id.csv'))
  
  # remain the minimun step of each from-to feature ##############################
  step_result <- do.call(rbind, step_result)
  step_result$rr_ID <- paste(step_result$From_feature,step_result$To_feature)
  
  step_result$rr_id_kegg <- sapply(seq(nrow(step_result)), function(i){
    rr_id <- c(step_result$from_id[i], 
               step_result$to_id[i])
    rr_id <- rr_id[order(rr_id)]
    rr_id <- paste0(rr_id, collapse = '_')
    return(rr_id)
  }, simplify = TRUE)
  
  # step_result_1 <- step_result[which(step_result$reaction_type == 'known'), ]
  # step_result_2 <- step_result[which(!step_result$reaction_type == 'known'), ]
  step_result_1 <- step_result
  rr_ID_unique <- unique(step_result_1$rr_ID)
  i <- 1
  step_result_i <- NULL
  step_result_i_min <- NULL
  step_result_min <- NULL
  for (i in 1:length(rr_ID_unique)) {
    step_result_i <-  step_result_1[ which(step_result_1$rr_ID==rr_ID_unique[i]),]
    
    if(length(which.min(step_result_i$step)) == 0){
      step_result_i_min <- step_result_i
    }else{
      step_result_i_min <- step_result_i[which.min(step_result_i$step),]
    }
    
    step_result_min <- rbind(step_result_min,step_result_i_min )
  }
  # step_result_min <- rbind(step_result_min, step_result_2)
  step_result_min$from_to_feature <- paste0(step_result_min$From_feature,
                                            '_', 
                                            step_result_min$To_feature)
  
  all_unique_from_to <- unique(step_result_min$from_to_feature)
  
  step_result_min <- lapply(all_unique_from_to, function(from_to){
    res <- step_result_min[which(step_result_min$from_to_feature == from_to),]
    res_combiend <- data.frame(From_feature = res$From_feature[1], 
                               To_feature = res$To_feature[1],
                               from_id = paste0(res$from_id, collapse = ';'), 
                               to_id = paste0(res$to_id, collapse = ';'), 
                               step = paste0(res$step, collapse = ';'), 
                               from_name = paste0(res$from_name, collapse = ';'), 
                               to_name = paste0(res$to_name, collapse = ';'), 
                               from_formula = paste0(res$from_formula, collapse = ';'), 
                               to_formula = paste0(res$to_formula, collapse = ';'), 
                               reaction_type = paste0(res$reaction_type, collapse = ';'))
    return(res_combiend)
  })
  step_result_min <- do.call(rbind, step_result_min)
  step_result_min$from_to_feature <- paste0(step_result_min$From_feature,
                                            '_', 
                                            step_result_min$To_feature)
  edge_table_add_id_min <- dplyr::left_join(edge_table, 
                                            step_result_min, 
                                            by = 'from_to_feature')
  edge_table_add_id_min$from_to_feature <- NULL
  edge_table_add_id_min$From_feature.y <- NULL
  edge_table_add_id_min$To_feature.y <- NULL
  # edge_table_add_id_min$rr_ID <- NULL
  # edge_table_add_id_min$rr_id_kegg <- NULL
  # edge_table_add_id_min$from_seq <- NULL
  # edge_table_add_id_min$to_seq <- NULL
  # edge_table_add_id_min$reaction_type[which(is.na(edge_table_add_id_min$reaction_type))] <- 'unknown_3'
  edge_table_add_id_min$reaction_type[which(is.na(edge_table_add_id_min$reaction_type))] <- 'Unknown2'
  
  colnames(edge_table_add_id_min)[which(colnames(edge_table_add_id_min) == 'From_feature.x')] <- 'From_feature'
  colnames(edge_table_add_id_min)[which(colnames(edge_table_add_id_min) == 'To_feature.x')] <- 'To_feature'
  idx_3 <- which(edge_table_add_id_min$From_feature %in% all_feature_in_pool & !edge_table_add_id_min$To_feature %in% all_feature_in_pool)
  idx_4 <- which(!edge_table_add_id_min$From_feature %in% all_feature_in_pool & edge_table_add_id_min$To_feature %in% all_feature_in_pool)
  
  
  idx_3_3 <- match(edge_table_add_id_min$From_feature[idx_3], id_pool$name)
  edge_table_add_id_min$from_id[idx_3] <- id_pool$all_id[idx_3_3]
  edge_table_add_id_min$from_name[idx_3] <- id_pool$all_name[idx_3_3]
  edge_table_add_id_min$from_formula[idx_3] <- id_pool$all_formula[idx_3_3]
  
  idx_4_4 <- match(edge_table_add_id_min$To_feature[idx_4], id_pool$name)
  edge_table_add_id_min$to_id[idx_4] <- id_pool$all_id[idx_4_4]
  edge_table_add_id_min$to_name[idx_4] <- id_pool$all_name[idx_4_4]
  edge_table_add_id_min$to_formula[idx_4] <- id_pool$all_formula[idx_4_4]
  
  all_idx_2 <- c(idx_3, 
               idx_4)
  
  # edge_table_add_id_min$reaction_type[all_idx_2] <- 'unknown_2'
  edge_table_add_id_min$reaction_type[all_idx_2] <- 'Unknown1'
  
  edge_table_add_id_min_final <- edge_table_add_id_min[, c('From', 'To', 'Type', 'MIDSimilarityScore', 'MotifSource', 'BestMotif', 'Reference', 'ReferenceRange', 
                                                           'ms2_score', 'From_feature', 'To_feature', 'DetlaMassAbs', 'DetlaMassAccuracy', 'AtomDiff', 
                                                           'ReactionClass', 'ReactionClassAll', 'ReactionFreq', 'step', 'reaction_type')]
  
  colnames(edge_table_add_id_min_final) <- c('From', 'To', 'Type', 'ISOSimilarity', 'MotifSource', 'BestMotif', 'Reference', 'ReferenceRange', 
                                             'MS2Similarity', 'From_feature', 'To_feature', 'DetlaMassAbs', 'DetlaMassAccuracy', 'AtomDiff', 
                                             'ReactionClass', 'ReactionClassAll', 'ReactionFreq', 'KEGGReactionStep', 'ReactionType')
  edge_table_add_id_min_final$From <- edge_table_add_id_min_final$From_feature
  edge_table_add_id_min_final$To <- edge_table_add_id_min_final$To_feature
  edge_table_add_id_min_final$From_feature <- NULL
  edge_table_add_id_min_final$To_feature <- NULL
  
  # modify the output path of final result
  dir.create(file.path(dir_path, 'Result'), showWarnings = FALSE, recursive = TRUE)
  
  # readr::write_csv(edge_table_add_id_min_final,
  #                  file = file.path(dir_path, 'Result','Labeled metabolite pairs.csv'))
  
  readr::write_csv(edge_table_add_id_min_final,
                   file = file.path(dir_path, '00_intermidate_data', 'Labeled metabolite pairs without RT filter.csv'))
  
  #### finally filter the edge annotation with RT ####
  from_rt_idx <- match(edge_table_add_id_min_final$From, id_table$name)
  from_rt <- as.numeric(id_table$rt[from_rt_idx])

  to_rt_idx <- match(edge_table_add_id_min_final$To, id_table$name)
  to_rt <- as.numeric(id_table$rt[to_rt_idx])
  
  rt_diff <- abs(from_rt - to_rt)
  rt_diff_idx <- which(rt_diff > 3)
  
  edge_table_add_id_min_final$DetlaMassAccuracy[-rt_diff_idx] <- NA
  edge_table_add_id_min_final$AtomDiff[-rt_diff_idx] <- NA
  edge_table_add_id_min_final$ReactionClass[-rt_diff_idx] <- NA
  edge_table_add_id_min_final$ReactionClassAll[-rt_diff_idx] <- NA
  edge_table_add_id_min_final$ReactionFreq[-rt_diff_idx] <- NA
  # edge_table_add_id_min_final$KEGGReactionStep[-rt_diff_idx] <- NA
  
  readr::write_csv(edge_table_add_id_min_final,
                   file = file.path(dir_path, 'Result','Labeled metabolite pairs.csv'))
}
