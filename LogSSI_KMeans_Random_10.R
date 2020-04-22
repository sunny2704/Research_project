library(stringr)
library(cluster)
library(proxy)
library(dtw)
library(dtwclust)
library(spatstat)

bf_c <- as.data.frame(read.csv("Wine_ClusterInfo.csv", header = F))
b_tr <- as.data.frame(read.csv("Wine_TRAIN.csv", header = F))
b_te <- as.data.frame(read.csv("Wine_TEST.csv", header = F))
full_data <- rbind.data.frame(b_tr, b_te)
full_dba <- as.numeric(dba(full_data, trace = FALSE))
grp_t <- 'INFO Group Value: '
grp_k_random <- 'INFO Cluster Method: Kmeans.Random'
grp_start <- 'INFO Start Random Kmeans Round'
grp_end <- 'INFO End  Random Kmeans Round: 10'
grp_dba <- "INFO Group DBA Value: "
k_dba <- ".*INFO Group DBA Value: "
k_ge <- ".*(INFO Group Elements: )"
kw <- ".*(INFO Group Value: )"
date_chk <- "2020-03-22"
loc_grp_v <- grep(grp_t, bf_c$V1, ignore.case = T)
BGSS <- 0.00
WGSS <- 0.00

for (i in 1:(length(loc_grp_v)))
{
  
  browser()
  if (i != (length(loc_grp_v)))
  {
    
    p <- q <- 1
    bcd <- list()
    act_ddc <- as.data.frame(bf_c$V1[loc_grp_v[i]: (loc_grp_v[(i+1)]-1)])
    colnames(act_ddc) <- "Values"
    start_pt <- grep(grp_start, act_ddc$Values, ignore.case = T)
    end_pt <- grep(grp_end, act_ddc$Values, ignore.case = T)
    
    
    for (z in 1:(length(start_pt)))
    {
      if(z != (length(start_pt)))
      {
        new_ddc <- as.data.frame(act_ddc$Values[start_pt[z]: (start_pt[(z+1)]-1)])
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
            val1 <- as.numeric(ddc_cls_dba)
            val2 <- as.numeric(full_dba)
            dtw_value_calc <- dtw(val1, val2,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
            BGSS <- BGSS + (y* (abs(dtw_value_calc$distance^2)))
            
          }
          if (length(ddc_clust) == 1)
          {
            new_1 <-  sub(k_ge, "", ddc_clust[1])
            bcd[y] <- list(new_1)
            
          }
        }
        
        p <- 1
        q <- 1
        
        for (p in 1: length(bcd))
        { 
          for (q in 1: length(bcd[[p]]))
          {
            u_val <- as.numeric(trimws(bcd[[p]][[q]]))
            v_val <- as.numeric(full_dba)
            d_n_val <- as.numeric(full_data[u_val,])
            dy_dist <- dtw(d_n_val,v_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
            WGSS <- WGSS + (dy_dist$distance^2)
          }
          
        }
        
        N <- 40
        
        answer <- log(BGSS/WGSS)
        #((N-K)/(K-1))* (BGSS/WGSS)
        
        lss <- round(answer,2)
        grp_val <- sub(kw, "", bf_c$V1[loc_grp_v[[i]]])
        grp_val <- str_trim(grp_val, "right")
        print(paste0("Group Number ", grp_val))
        print(paste0("Round Number ", z))
        print("Clustering Type: Random")
        print(paste0("The Log SS Ratio index is ", lss))
        BGSS <- 0.00
        WGSS <- 0.00
        
      }
      if(z == (length(start_pt)))
      {
        
        new_ddc <- as.data.frame(act_ddc$Values[start_pt[z]: nrow(act_ddc)])
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
            val1 <- as.numeric(ddc_cls_dba)
            val2 <- as.numeric(full_dba)
            dtw_value_calc <- dtw(val1, val2,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
            BGSS <- BGSS + (y* (abs(dtw_value_calc$distance^2)))
            
          }
          if (length(ddc_clust) == 1)
          {
            new_1 <-  sub(k_ge, "", ddc_clust[1])
            bcd[y] <- list(new_1)
            
          }
        }
        
        p <- 1
        q <- 1
        
        for (p in 1: length(bcd))
        { 
          for (q in 1: length(bcd[[p]]))
          {
            u_val <- as.numeric(trimws(bcd[[p]][[q]]))
            v_val <- as.numeric(full_dba)
            d_n_val <- as.numeric(full_data[u_val,])
            dy_dist <- dtw(d_n_val,v_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
            WGSS <- WGSS + (dy_dist$distance^2)
          }
          
        }
        
        N <- 40
        
        answer <- log(BGSS/WGSS)
        #((N-K)/(K-1))* (BGSS/WGSS)
        
        lss <- round(answer,2)
        grp_val <- sub(kw, "", bf_c$V1[loc_grp_v[[i]]])
        grp_val <- str_trim(grp_val, "right")
        print(paste0("Group Number ", grp_val))
        print(paste0("Round Number ", z))
        print("Clustering Type: Random")
        print(paste0("The Log SS Ratio index is ", lss))
        BGSS <- 0.00
        WGSS <- 0.00
        
        
      }
    }
    
  }
  
  
  if (i == length(loc_grp_v))
  {
    p <- q <- 1
    BGSS <- 0.00
    WGSS <- 0.00
    bcd <- list()
    
    act_ddc <- as.data.frame(bf_c$V1[loc_grp_v[[i]]:nrow(bf_c)])
    colnames(act_ddc) <- "Values"
    
    start_pt <- grep(grp_start, act_ddc$Values, ignore.case = T)
    end_pt <- grep(grp_end, act_ddc$Values, ignore.case = T)
    
    BGSS <- 0.00
    WGSS <- 0.00
    
    for (z in 1:(length(start_pt)))
    {
      if(z != (length(start_pt)))
      {
        new_ddc <- as.data.frame(act_ddc$Values[start_pt[z]: (start_pt[(z+1)]-1)])
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
            
            val1 <- as.numeric(ddc_cls_dba)
            val2 <- as.numeric(full_dba)
            dtw_value_calc <- dtw(val1, val2,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
            BGSS <- BGSS + (y* (abs(dtw_value_calc$distance^2)))
            
          }
          if (length(ddc_clust) == 1)
          {
            new_1 <-  sub(k_ge, "", ddc_clust[1])
            bcd[y] <- list(new_1)
            
          }
        }   
        
        
        p <- 1
        q <- 1
        
        for (p in 1: length(bcd))
        { 
          for (q in 1: length(bcd[[p]]))
          {
            u_val <- as.numeric(trimws(bcd[[p]][[q]]))
            v_val <- as.numeric(full_dba)
            d_n_val <- as.numeric(full_data[u_val,])
            dy_dist <- dtw(d_n_val,v_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
            WGSS <- WGSS + (dy_dist$distance^2)
          }
          
        }
        
        N <- 40
        
        answer <- log(BGSS/WGSS)
        #((N-K)/(K-1))* (BGSS/WGSS)
        
        lss <- round(answer,2)
        grp_val <- sub(kw, "", bf_c$V1[loc_grp_v[[i]]])
        grp_val <- str_trim(grp_val, "right")
        print(paste0("Group Number ", grp_val))
        print(paste0("Round Number ", z))
        print("Clustering Type: Random")
        print(paste0("The Log SS Ratio index is ", lss))
        BGSS <- 0.00
        WGSS <- 0.00
        
      }
      if(z == (length(start_pt)))
      {
        
        new_ddc <- as.data.frame(act_ddc$Values[start_pt[z]: nrow(act_ddc)])
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
            
            val1 <- as.numeric(ddc_cls_dba)
            val2 <- as.numeric(full_dba)
            dtw_value_calc <- dtw(val1, val2,keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
            BGSS <- BGSS + (y* (abs(dtw_value_calc$distance^2)))
            
          }
          if (length(ddc_clust) == 1)
          {
            new_1 <-  sub(k_ge, "", ddc_clust[1])
            bcd[y] <- list(new_1)
            
          }
        }   
        
        
        p <- 1
        q <- 1
        
        for (p in 1: length(bcd))
        { 
          for (q in 1: length(bcd[[p]]))
          {
            u_val <- as.numeric(trimws(bcd[[p]][[q]]))
            v_val <- as.numeric(full_dba)
            d_n_val <- as.numeric(full_data[u_val,])
            dy_dist <- dtw(d_n_val,v_val, keep=TRUE, step.pattern = symmetric1, open.end = FALSE)
            WGSS <- WGSS + (dy_dist$distance^2)
          }
          
        }
        
        N <- 40
        
        answer <- log(BGSS/WGSS)
        #((N-K)/(K-1))* (BGSS/WGSS)
        
        lss <- round(answer,2)
        grp_val <- sub(kw, "", bf_c$V1[loc_grp_v[[i]]])
        grp_val <- str_trim(grp_val, "right")
        print(paste0("Group Number ", grp_val))
        print(paste0("Round Number ", z))
        print("Clustering Type: Random")
        print(paste0("The Log SS Ratio index is ", lss))
        BGSS <- 0.00
        WGSS <- 0.00
      }
    }
    
  }
}

