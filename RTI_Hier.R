library(stringr)
library(cluster)
library(proxy)
library(dtw)
library(dtwclust)
library(spatstat)
library(clusterCrit)

bf_c <- as.data.frame(read.csv("Gun_Point_ClusterInfo.csv", header = F))
b_tr <- as.data.frame(read.csv("Gun_Point_TRAIN.csv", header = F))
b_te <- as.data.frame(read.csv("Gun_Point_TEST.csv", header = F))
full_data <- rbind.data.frame(b_tr, b_te)
full_dba <- as.numeric(dba(full_data, trace = FALSE))
grp_t <- 'INFO Group Value: '
grp_k_random <- 'INFO Cluster Method: Hierarchy.Complete'
grp_dba <- "INFO Group DBA Value: "
k_ge <- ".*(INFO Group Elements: )"
kw <- ".*(INFO Group Value: )"
k_dba <- ".*INFO Group DBA Value: "
date_chk <- "2020-03-22"
loc_grp_v <- grep(grp_t, bf_c$V1, ignore.case = T)
# BGSS <- 0.00
# WGSS <- 0.00

for (i in 1:(length(loc_grp_v)))
{
  if (i != (length(loc_grp_v)))
  {
    WGSS <- 0.00
    min_db <- 0.00
    p <- q <- 1
    bcd <- list()
    bcd_dba <- list()
    new_ddc <- as.data.frame(bf_c$V1[loc_grp_v[i]: (loc_grp_v[(i+1)]-1)])
    colnames(new_ddc) <- "Values"
    loc_ddc <- grep(grp_k_random, new_ddc$Values, ignore.case = T)
    K <- length(loc_ddc)
    
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
      
      new_f_dba = sub(k_dba, "", dba_data[1])
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
        new_ddc <- as.data.frame(new_ddc)
        colnames(new_ddc) <- "Values"
        bcd[y] <- list(ddc_cluster)
        # val1 <- as.numeric(ddc_cls_dba)
        # val2 <- as.numeric(full_dba)
        # dtw_value_calc <- dtw(val1, val2,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        # BGSS <- BGSS + (y* (abs(dtw_value_calc$distance)^2))
        
      }
      if (length(ddc_clust) == 1)
      {
        ddc_cluster <-  sub(k_ge, "", ddc_clust[1])
        bcd[y] <- list(ddc_cluster)
        
      }
      
      for(x in 1:length(ddc_cluster))
      {
        val_a <- as.numeric(ddc_cls_dba)
        val_b <- as.numeric(trimws(ddc_cluster[x]))
        b_val <-as.numeric(full_data[val_b,])
        dtw_loc_dba <- dtw(val_a, b_val,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        WGSS <- WGSS + (abs(dtw_loc_dba$distance)^2)
        
      }
    }
    
    p <- 1
    q <- 1
    
    min_db = 99999999
    
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
          chk_min_dist <- (abs(dy_dist$distance)^2)
          if(chk_min_dist < min_db)
          {
            min_db <- chk_min_dist
          }
        }
        
      }
      
    }
    
    
    
    
    N <- 40
    
    
    ans <- WGSS/min_db
    answer <- ans/N
    
    #((N-K)/(K-1))* (BGSS/WGSS)
    
    rti <- round(answer,2)
    grp_val <- sub(kw, "", bf_c$V1[loc_grp_v[[i]]])
    grp_val <- str_trim(grp_val, "right")
    print(paste0("Group Number ", grp_val))
    print("Clustering Type: Hierarchy")
    print(paste0("The Ray Turi index is ", rti))
    WGSS <- 0.00
    min_db <- 0.00
    
  }
  
  
  if (i == length(loc_grp_v))
  {
    p <- q <- 1
    WGSS <- 0.00
    min_db <- 0.00
    bcd <- list()
    bcd_dba <- list()
    
    new_ddc <- as.data.frame(bf_c$V1[loc_grp_v[[i]]:nrow(bf_c)])
    colnames(new_ddc) <- "Values"
    loc_ddc <- grep(grp_k_random, new_ddc$Values, ignore.case = T)
    K <- length(loc_ddc)
    
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
      
      new_f_dba = sub(k_dba, "", dba_data[1])
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
        new_ddc <- as.data.frame(new_ddc)
        colnames(new_ddc) <- "Values"
        bcd[y] <- list(ddc_cluster)
        # 
        # val1 <- as.numeric(ddc_cls_dba)
        # val2 <- as.numeric(full_dba)
        # dtw_value_calc <- dtw(val1, val2,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        # BGSS <- BGSS + (y* (abs(dtw_value_calc$distance)^2))
        
      }
      if (length(ddc_clust) == 1)
      {
        ddc_cluster <-  sub(k_ge, "", ddc_clust[1])
        bcd[y] <- list(new_1)
        
      }
      
      for(x in 1:length(ddc_cluster))
      {
        val_a <- as.numeric(ddc_cls_dba)
        val_b <- as.numeric(trimws(ddc_cluster[x]))
        b_val <-as.numeric(full_data[val_b,])
        dtw_loc_dba <- dtw(val_a, b_val,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        WGSS <- WGSS + (abs(dtw_loc_dba$distance)^2)
        
      }
      
    }   
    
    
    p <- 1
    q <- 1
    min_db = 99999999
    
    
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
          chk_min_dist <- (abs(dy_dist$distance)^2)
          if(chk_min_dist < min_db)
          {
            min_db <- chk_min_dist
          }
        }
        
      }
      
    }
    
    
    
    
    
    
    N <- 40
    
    ans <- WGSS/min_db
    answer <- ans/N
    #((N-K)/(K-1))* (BGSS/WGSS)
    
    rti <- round(answer,2)
    grp_val <- sub(kw, "", bf_c$V1[loc_grp_v[[i]]])
    grp_val <- str_trim(grp_val, "right")
    print(paste0("Group Number ", grp_val))
    print("Clustering Type: Hierarchy")
    print(paste0("The Ray Turi index is ", rti))
    WGSS <- 0.00
    min_db <- 0.00
    
  }
}

