library(stringr)
library(cluster)
library(proxy)
library(dtw)
library(dtwclust)
library(spatstat)

bf_c <- as.data.frame(read.csv("ToeSegmentation2_ClusterInfo.csv", header = F))
b_tr <- as.data.frame(read.csv("ToeSegmentation2_TRAIN.csv", header = F))
b_te <- as.data.frame(read.csv("ToeSegmentation2_TEST.csv", header = F))
full_data <- rbind.data.frame(b_tr, b_te)
full_dba <- as.numeric(dba(full_data, trace = FALSE))
grp_t <- 'INFO Group Value: '
# grp_ddc <- 'INFO Cluster Method: DDC'
grp_k_random <- 'INFO Cluster Method: Hierarchy.Single'
grp_dba <- "INFO Group DBA Value: "
k_ge <- ".*(INFO Group Elements: )"
kw <- ".*INFO Group DBA Value: "
k_new <- ".*(INFO Group Value: )"
date_chk <- "2020-03-22"
loc_grp_v <- grep(grp_t, bf_c$V1, ignore.case = T)
# BGSS <- 0.00
# WGSS <- 0.00

for (i in 1:(length(loc_grp_v)))
{
  if (i != (length(loc_grp_v)))
  {
    
    p <- q <- 1
    bcd <- list()
    bcd_dba <- list()
    new_ddc <- as.data.frame(bf_c$V1[loc_grp_v[i]: (loc_grp_v[(i+1)]-1)])
    colnames(new_ddc) <- "Values"
    loc_ddc <- grep(grp_k_random, new_ddc$Values, ignore.case = T)
    K <- length(loc_ddc)
    
    
    sum_loc_dba <- 0.00
    sum_full_dba <- 0.00
    for (y in 1:length(loc_ddc))
    {
      end <- nrow(new_ddc)
      new_ddc_sample <-as.character(new_ddc$Values[(loc_ddc[y]+2) :end])
      
      dba_set <- as.data.frame(new_ddc$Values[(loc_ddc[y]+3):end])
      colnames(dba_set) <- "Values"
      dba_loc <- grep(grp_dba, dba_set$Values, ignore.case = T)
      date_loc <- grep(date_chk, dba_set$Values, ignore.case = T)
      #local_dba
      dba_data <- as.character(dba_set$Values[(date_loc[[1]]): (date_loc[[2]] - 1)])
      
      new_f_dba = sub(kw, "", dba_data[1])
      ddc_dba_len <- length(dba_data)
      ddc_dba1 = dba_data[2:ddc_dba_len]
      ddc_cls_dba <- c(new_f_dba, ddc_dba1)
      bcd_dba[y] <- list(ddc_cls_dba)
      
      loc_dba <- grep(grp_dba, new_ddc_sample)
      ddc_clust <- new_ddc_sample[1:(loc_dba[[1]] - 1)]
      if (length(ddc_clust) != 1)
      {
        new_1 <-  sub(k_ge, "", ddc_clust[1])
        ddc_clust_len <- length(ddc_clust)
        ddc_clust1 = ddc_clust[2:ddc_clust_len]
        ddc_cluster <- c(new_1, ddc_clust1)

        
      }
      if (length(ddc_clust) == 1)
      {
        ddc_cluster <-  sub(k_ge, "", ddc_clust[1])
        # bcd[y] <- list(ddc_cluster)
        
      }
      
      for(x in 1:length(ddc_cluster))
      {
        val_a <- as.numeric(ddc_cls_dba)
        val_b <- as.numeric(full_dba)
        val_c <- as.numeric(trimws(ddc_cluster[x]))
        c_val <-as.numeric(full_data[val_c,])
        dtw_loc_dba <- dtw(val_a, c_val,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        dtw_full_dba <- dtw(val_b, c_val,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        sum_loc_dba <- sum_loc_dba + dtw_loc_dba$distance
        sum_full_dba <- sum_full_dba + dtw_full_dba$distance
      }
      
    }
    
    p <- 1
    q <- 1
    max_db = 0.00
    
    for (p in 1: length(bcd_dba))
    { 
      # print(p)
      for (q in 1: length(bcd_dba))
      {
        # print(q)
        if(p!= q)
        {
          u_val <- as.numeric(unlist(bcd_dba[p]))
          v_val <- as.numeric(unlist(bcd_dba[q]))
          # d_n_val <- as.numeric(full_data[u_val,])
          dy_dist <- dtw(u_val,v_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
          if(dy_dist$distance > max_db)
          {
            max_db <- dy_dist$distance
          }
        }
        
      }
      
    }
    
    parta <- max_db * sum_full_dba
    partb <- K * sum_loc_dba
    ans <- parta/ partb
    
    # ans <- (max_db/K)*(sum_full_dba/sum_loc_dba)
    answer <- ans^2
    
    
    pbm <- round(answer,2)
    grp_val <- sub(k_new, "", bf_c$V1[loc_grp_v[[i]]])
    grp_val <- str_trim(grp_val, "right")
    print(paste0("Group Number ", grp_val))
    print("Clustering Type: Hierarchy")
    print(paste0("The PBM index is ", pbm))
    sum_loc_dba <- 0.00
    sum_full_dba <- 0.00
    max_db <- 0.00
    
  }
  
  
  if (i == length(loc_grp_v))
  {
    p <- q <- 1
    bcd <- list()
    bcd_dba <- list()
    
    new_ddc <- as.data.frame(bf_c$V1[loc_grp_v[[i]]:nrow(bf_c)])
    colnames(new_ddc) <- "Values"
    loc_ddc <- grep(grp_k_random, new_ddc$Values, ignore.case = T)
    K <- length(loc_ddc)
    sum_loc_dba <- 0.00
    sum_full_dba <- 0.00
    

    
    for (y in 1:length(loc_ddc))
    { 
      
      end <- nrow(new_ddc)
      new_ddc_sample <-as.character(new_ddc$Values[(loc_ddc[y]+2) :end])
      
      dba_set <- as.data.frame(new_ddc$Values[(loc_ddc[y]+3):end])
      colnames(dba_set) <- "Values"
      dba_loc <- grep(grp_dba, dba_set$Values, ignore.case = T)
      date_loc <- grep(date_chk, dba_set$Values, ignore.case = T)
      #local_dba
      dba_data <- as.character(dba_set$Values[(date_loc[[1]]): (date_loc[[2]] - 1)])
      
      new_f_dba = sub(kw, "", dba_data[1])
      ddc_dba_len <- length(dba_data)
      ddc_dba1 = dba_data[2:ddc_dba_len]
      ddc_cls_dba <- c(new_f_dba, ddc_dba1)
      bcd_dba[y] <- list(ddc_cls_dba)
      
      loc_dba <- grep(grp_dba, new_ddc_sample)
      ddc_clust <- new_ddc_sample[1:(loc_dba[[1]] - 1)]
      
      if (length(ddc_clust) != 1)
      {
        new_1 <-  sub(k_ge, "", ddc_clust[1])
        ddc_clust_len <- length(ddc_clust)
        ddc_clust1 = ddc_clust[2:ddc_clust_len]
        ddc_cluster <- c(new_1, ddc_clust1)
        
      }
      if (length(ddc_clust) == 1)
      {
        ddc_cluster  <-  sub(k_ge, "", ddc_clust[1])
        bcd[y] <- list(ddc_cluster)
        
      }
      
      for(x in 1:length(ddc_cluster))
      {
        val_a <- as.numeric(ddc_cls_dba)
        val_b <- as.numeric(full_dba)
        val_c <- as.numeric(trimws(ddc_cluster[x]))
        c_val <-as.numeric(full_data[val_c,])
        dtw_loc_dba <- dtw(val_a, c_val,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        dtw_full_dba <- dtw(val_b, c_val,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        sum_loc_dba <- sum_loc_dba + dtw_loc_dba$distance
        sum_full_dba <- sum_full_dba + dtw_full_dba$distance
      }
      
    }   
    
    
    p <- 1
    q <- 1
    max_db = 0.00
    
    for (p in 1: length(bcd_dba))
    { 
      for (q in 1: length(bcd_dba))
      {
        if(p!= q)
        {
          u_val <- as.numeric(unlist(bcd_dba[p]))
          v_val <- as.numeric(unlist(bcd_dba[q]))
          # d_n_val <- as.numeric(full_data[u_val,])
          dy_dist <- dtw(u_val,v_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
          if(dy_dist$distance > max_db)
          {
            max_db <- dy_dist$distance
          }
        }
        
      }
      
    }
    
    parta <- max_db * sum_full_dba
    partb <- K * sum_loc_dba
    ans <- parta/ partb
    # ans <- (max_db/K)*(sum_full_dba/sum_loc_dba)
    answer <- ans^2
    
    
    pbm <- round(answer,2)
    grp_val <- sub(k_new, "", bf_c$V1[loc_grp_v[[i]]])
    grp_val <- str_trim(grp_val, "right")
    print(paste0("Group Number ", grp_val))
    print("Clustering Type: Hierarchy")
    print(paste0("The PBM index is ", pbm))
    sum_loc_dba <- 0.00
    sum_full_dba <- 0.00
    max_db <- 0.00
    
  }
}

