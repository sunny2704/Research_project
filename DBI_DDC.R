library(stringr)
library(cluster)
library(proxy)
library(dtw)
library(dtwclust)
library(spatstat)
library(janitor)

bf_c <- as.data.frame(read.csv("ToeSegmentation2_ClusterInfo.csv", header = F))
b_tr <- as.data.frame(read.csv("ToeSegmentation2_TRAIN.csv", header = F))
b_te <- as.data.frame(read.csv("ToeSegmentation2_TEST.csv", header = F))
full_data <- rbind.data.frame(b_tr, b_te)
full_dba <- as.numeric(dba(full_data, trace = FALSE))
grp_t <- 'INFO Group Value: '
grp_ddc <- 'INFO Cluster Method: DDC'
grp_dba <- "INFO Group DBA Value: "
k_ge <- ".*(INFO Group Elements: )"
k_dba <- ".*(INFO Group DBA Value: )"
kw <- ".*(INFO Group Value: )"
date_chk <- "2019-02-17"
loc_grp_v <- grep(grp_t, bf_c$V1, ignore.case = T)
BGSS <- 0.00
WGSS <- 0.00

for (i in 1:(length(loc_grp_v)))
{
  if (i != (length(loc_grp_v)))
  {
    
    p <- q <- 1
    bcd <- list()
    dba_list <- list()
    new_ddc <- as.data.frame(bf_c$V1[loc_grp_v[i]: (loc_grp_v[(i+1)]-1)])
    colnames(new_ddc) <- "Values"
    loc_ddc <- grep(grp_ddc, new_ddc$Values, ignore.case = T)
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
      dba_list[y] <-list(ddc_cls_dba)
      
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
        
      }
      if (length(ddc_clust) == 1)
      {
        new_1 <-  sub(k_ge, "", ddc_clust[1])
        bcd[y] <- list(new_1)
        
      }
    }
    
    p <- 1
    q <- 1
    temp_str <- setNames(data.frame(matrix(ncol = 3, nrow = 200), stringsAsFactors = FALSE), c("cluster_no", "mean_dist","dba_val"))
    
    
    for (p in 1: length(bcd))
    { 
      for (q in 1: length(bcd[[p]]))
      {
        u_val <- as.numeric(trimws(bcd[[p]][[q]]))
        v_val <- as.numeric(trimws(dba_list[[p]]))
        d_n_val <- as.numeric(full_data[u_val,])
        dy_dist <- dtw(d_n_val,v_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        WGSS <- WGSS + dy_dist$distance
      }
      temp_str$cluster_no[p] <- p
      temp_str$mean_dist[p] <- WGSS/(length(bcd[[p]]))
      temp_var <- v_val
      temp_str$dba_val[p] <- list(temp_var)
      
    }
    
    
    
    temp_str <- remove_empty(temp_str,which = "rows")
    numerator <- 0.00
    for (f in 1:nrow(temp_str))
    {
      f_val <-as.numeric(unlist(temp_str$dba_val[f]))
      
      fg_dist <- 0.00
      temp_value <- 0.00
      for (g in 1:nrow(temp_str))
      {
        if(f !=g)
        {
          g_val <- as.numeric(unlist(temp_str$dba_val[g]))
          fg_dtw <- dtw(f_val,g_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
          temp_value <- (temp_str$mean_dist[p] + temp_str$mean_dist[p])/ fg_dtw$distance
          if(temp_value > fg_dist)
          {
            fg_dist <- temp_value
          }
        }
      }
    
      
      numerator <- numerator + temp_value
    }
    
    
    answer <- numerator / length(loc_ddc)
    #((N-K)/(K-1))* (BGSS/WGSS)
    
    dbi <- round(answer,2)
    grp_val <- sub(kw, "", bf_c$V1[loc_grp_v[[i]]])
    grp_val <- str_trim(grp_val, "right")
    print(paste0("Group Number ", grp_val))
    print("Clustering Type: DDC")
    print(paste0("The The Davies-Bouldin index is ", dbi))
    BGSS <- 0.00
    WGSS <- 0.00
    
  }
  
  
  if (i == length(loc_grp_v))
  {
    p <- q <- 1
    BGSS <- 0.00
    WGSS <- 0.00
    bcd <- list()
    
    new_ddc <- as.data.frame(bf_c$V1[loc_grp_v[[i]]:nrow(bf_c)])
    colnames(new_ddc) <- "Values"
    loc_ddc <- grep(grp_ddc, new_ddc$Values, ignore.case = T)
    K <- length(loc_ddc)
    BGSS <- 0.00
    WGSS <- 0.00
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
      dba_list[y] <-list(ddc_cls_dba)
      
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
        
      }
      if (length(ddc_clust) == 1)
      {
        new_1 <-  sub(k_ge, "", ddc_clust[1])
        bcd[y] <- list(new_1)
        
      }
    }   
    
    
    p <- 1
    q <- 1
    temp_str <- setNames(data.frame(matrix(ncol = 3, nrow = 200), stringsAsFactors = FALSE), c("cluster_no", "mean_dist","dba_val"))
    
    
    for (p in 1: length(bcd))
    { 
      for (q in 1: length(bcd[[p]]))
      {
        u_val <- as.numeric(trimws(bcd[[p]][[q]]))
        v_val <- as.numeric(trimws(dba_list[[p]]))
        d_n_val <- as.numeric(full_data[u_val,])
        dy_dist <- dtw(d_n_val,v_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
        WGSS <- WGSS + dy_dist$distance
      }
      temp_str$cluster_no[p] <- p
      temp_str$mean_dist[p] <- WGSS/(length(bcd[[p]]))
      temp_var <- v_val
      temp_str$dba_val[p] <- list(temp_var)
      
    }
    
    # done till here!!!!
    
    
    temp_str <- remove_empty(temp_str,which = "rows")
    numerator <- 0.00
    for (f in 1:length(temp_str))
    {
      
      f_val <-as.numeric(unlist(temp_str$dba_val[f]))
      
      fg_dist <- 0.00
      temp_value <- 0.00
      
      for (g in 1:length(temp_str))
      {
        if(f !=g)
        {
          g_val <- as.numeric(unlist(temp_str$dba_val[g]))
          fg_dtw <- dtw(f_val,g_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
          temp_value <- (temp_str$mean_dist[p] + temp_str$mean_dist[p])/ fg_dtw$distance
          if(temp_value > fg_dist)
          {
            fg_dist <- temp_value
          }
        }
      }
      
      
      numerator <- numerator + temp_value
    }
    
    
    answer <- numerator / length(loc_ddc)
    #((N-K)/(K-1))* (BGSS/WGSS)
    
    dbi <- round(answer,2)
    grp_val <- sub(kw, "", bf_c$V1[loc_grp_v[[i]]])
    grp_val <- str_trim(grp_val, "right")
    print(paste0("Group Number ", grp_val))
    print("Clustering Type: DDC")
    print(paste0("The Davies-Bouldin index is ", dbi))
    BGSS <- 0.00
    WGSS <- 0.00
    
  }
}

