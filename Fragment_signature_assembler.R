if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChemmineR")

pascalTriangle_1 <- function(h) {
  lapply(h, function(i) choose(i, 1:i))
}
install.packages("lsa")

library("ChemmineR")
library(lsa)
mat_list_da <-as.matrix(read.csv(file.choose(), header=T, sep=",")) #MSMS file position and assigned biotransformation reactions data
mat_list_da_num <- length(mat_list_da[,1])

ppm_value <- 5
parent_site_uni <- 1
int_cut_level <- 0

#input pd structure
pd_structure_smile <- read.SDFset("D:/Structure2D_CID_20393.sdf") 
plot(pd_structure_smile, atomnum = TRUE ,no_print_atoms = "",atomcex=0.7, print=FALSE)

mat_raw_charge <- "CHARGE=1-"

mat_raw_ms1_mz <- matrix(nrow = mat_list_da_num ,ncol = 1000)
mat_raw_ms1_int <- matrix(nrow = mat_list_da_num ,ncol = 1000)
mat_raw_ms2_mz <- matrix(nrow = mat_list_da_num ,ncol = 1000)
mat_raw_ms2_int <- matrix(nrow = mat_list_da_num ,ncol = 1000)

i=1
for(i in c(1:mat_list_da_num)){
  mat_pathway <- mat_list_da[i,2]
  mat_list_da_n <- mat_list_da[i,1]
  mat_raw <- readLines(mat_pathway) 
  mat_raw_num <- length(mat_raw)
  j=1
  mat_rawmz1 <- unlist(strsplit(as.character(mat_raw[3]),": "))
  mat_rawmz <- mat_rawmz1[2]
  if(i==1){
    mat_rawmz_pd <- as.numeric(mat_rawmz)
  }else{
    
  }
  record_ms1 <- 0
  record_ms2 <- 0
  cord_ms1 <- 1
  cord_ms2 <- 1
  for(j in c(1:mat_raw_num)){
    
    if(mat_raw[j] == "IONTYPE: Positive\t"){
      mat_raw_charge <- "CHARGE=1+"
    }else if(mat_raw[j] == "IONTYPE: Negative\t"){
      mat_raw_charge <- "CHARGE=1-"
    }else{
      
    }
    mat_list_damsc <- unlist(strsplit(as.character(mat_raw[j]),"\t"))
    mat_list_damsc1 <- unlist(strsplit(as.character(mat_raw[j-2]),"\t"))
    
    
    if(j >= 3 && mat_list_damsc1[1] == "MSTYPE: MS1"){
      record_ms1 <- 1
    }else if(j >= 3 && mat_list_damsc[1] == "MSTYPE: MS2"){
      record_ms1 <- 0
    }else if(j >= 3 && mat_list_damsc1[1] == "MSTYPE: MS2"){
      record_ms2 <- 1
    }else{
      
    }
    
    if(record_ms1 == 1){
      mat_raw_ms1_rawvalue <- unlist(strsplit(as.character(mat_raw[j]),"\t"))
      mat_raw_ms1_mz[i,cord_ms1] <- mat_raw_ms1_rawvalue[1] 
      mat_raw_ms1_int[i,cord_ms1] <- mat_raw_ms1_rawvalue[2]
      cord_ms1 = cord_ms1 +1
    }else if(record_ms2 == 1){
      mat_raw_ms2_rawvalue <- unlist(strsplit(as.character(mat_raw[j]),"\t"))
      mat_raw_ms2_mz[i,cord_ms2] <- mat_raw_ms2_rawvalue[1] 
      mat_raw_ms2_int[i,cord_ms2] <- mat_raw_ms2_rawvalue[2]
      cord_ms2 = cord_ms2 +1
    }else{
      
    }
  }
}
pd_ms2_num <- length(na.omit(mat_raw_ms2_mz[1,]))

#assign structure of fragment
mat_raw_ms2_pd_structure <- matrix(nrow = pd_ms2_num ,ncol = 100)
mat_raw_ms2_pd_structure[,1] <- na.omit(mat_raw_ms2_mz[1,])

#manual input identified fragments of parent compound
mat_raw_ms2_pd_structure[1,2:11] <-c(3,4,14,15,16,17,18,19,20,42)
mat_raw_ms2_pd_structure[2,2:10] <-c(1,5,6,7,8,9,10,11,12)
mat_raw_ms2_pd_structure[3,2:11] <-c(1,2,9,13,14,15,16,17,18,19)
mat_raw_ms2_pd_structure[4,2:13] <-c(2,3,4,13,14,15,16,17,18,19,20,42)
mat_raw_ms2_pd_structure[5,2:16] <-c(1,2,5,6,8,9,10,12,13,14,15,16,17,18,19)
mat_raw_ms2_pd_structure[6,2:18] <-c(1,2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
mat_raw_ms2_pd_structure[7,2:22] <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,42)
fragment_list <- c(1,2,3,4,5,6,7)
fragment_mz_list <- matrix(nrow = length(fragment_list) ,ncol = 1)
i=1
for(i in c(1:length(fragment_list))){
  fragment_mz_list[i,1] <- mat_raw_ms2_mz[1,fragment_list[i]]
}

mat_raw_ms2_pd_structure <- mat_raw_ms2_pd_structure[complete.cases(mat_raw_ms2_pd_structure[,2]),]

#start calculate 

metabolite_id_result_table_score <- matrix(nrow = mat_list_da_num ,ncol = 10)
metabolite_id_result_table_score_p <- matrix(nrow = mat_list_da_num ,ncol = 10)

metabolite_id_result_table_score_p_1 <- matrix(nrow = mat_list_da_num ,ncol = 10)

cc=1
for(cc in c(1:10)){
  metabolite_id_result_posit_name <-  paste0("metabolite_id_result_posit_f_",cc)
  metabolite_id_result_table_name <-  paste0("metabolite_id_result_table_f_",cc)
  
  mat_raw_ms2_mz_pd <- mat_raw_ms2_mz[1,]
  mat_raw_ms2_mz_pd <- na.omit(mat_raw_ms2_mz_pd)
  mat_raw_ms2_mz_pd_id <- matrix(nrow = length(fragment_list) ,ncol = 1)
  metabolite_id_result_table <- matrix(nrow = length(fragment_list)+1 ,ncol = 1)
  metabolite_id_result_table[1,1] <- "name"
  df=1
  metabolite_id_result_posit <- matrix(nrow = 2 ,ncol = 1)
  metabolite_id_result_posit[1,1] <- "name"
  metabolite_id_result_posit[2,1] <- "site"
  i=1
  for(i in c(1:length(fragment_list))){
    metabolite_id_result_table[i+1,1] <- fragment_mz_list[i,1]
    mat_raw_ms2_mz_pd_id[i,1] <- as.numeric(mat_raw_ms2_mz_pd[fragment_list[i]])
  }
  
  reaction_type_table <- matrix(0, nrow = (mat_list_da_num-1) ,ncol = 5)
  i=2
  for(i in c(2:mat_list_da_num)){
    mat_pd_structure_ass <- matrix(nrow = length(fragment_list) ,ncol = 100)
    mat_pd_structure_ass[1:length(fragment_list),1] <- mat_raw_ms2_mz_pd_id[1:length(fragment_list),1]
    reaction_type <- unlist(strsplit(as.character(mat_list_da[i,(cc+2)]),"+", fixed = TRUE))
    reaction_type_table[i-1,1:length(reaction_type)] <- reaction_type
    reaction_num_total <- sum(as.numeric(unlist(pascalTriangle_1(length(reaction_type)))))
    reaction_type_table_1 <- matrix(nrow = reaction_num_total ,ncol = 5)
    s=1
    cc_count_x <- 1
    for(s in c(1:length(reaction_type))){
      kkk <- unique(t(combn(reaction_type,s)))
      j=1
      for(j in c(1:nrow(kkk))){
        reaction_type_table_1[cc_count_x, 1:length(kkk[j,])]<- kkk[j,]
        cc_count_x <- cc_count_x + 1
      }
    }
    reaction_type_table_1 <- reaction_type_table_1[complete.cases(reaction_type_table_1[,1]),]
    reaction_num_1 <- nrow(reaction_type_table_1)
    
    if(is.null(reaction_num_1)){
      reaction_num_1 = 1
    }else{
      
    }
    
    dd=1
    ds=1
    for(ds in c(1:reaction_num_1)){
      if(reaction_num_1 == 1){
        reaction_posit <- which(meta_comb_list_f[,2] == paste(na.omit(reaction_type_table_1[ds]),collapse="+"))
      }else{
        reaction_posit <- which(meta_comb_list_f[,2] == paste(na.omit(reaction_type_table_1[ds,]),collapse="+"))
      }
      reaction_mz_vzlue <- round(as.numeric(meta_comb_list_f[reaction_posit,3]),4)-as.numeric(mat_raw_ms1_mz[1,1])
      for(dd in c(1:length(fragment_list))){
        mat_pd_structure_ass[dd,ds+1] <- round(as.numeric(mat_pd_structure_ass[dd,1]),4)+ reaction_mz_vzlue
      }
    }
    mat_pd_structure_ass_1 <- t(mat_pd_structure_ass)
    mat_pd_structure_ass_1 <- mat_pd_structure_ass_1[complete.cases(mat_pd_structure_ass_1[,1]),]
    if(is.null(nrow(mat_pd_structure_ass_1))){
      mat_pd_structure_ass_1 <- matrix(unlist(mat_pd_structure_ass_1), ncol = length(fragment_list), byrow = TRUE)
    }else{
      
    }
    mat_pd_structure_ass <- t(mat_pd_structure_ass_1)
    rm(mat_pd_structure_ass_1)
    
    mat_raw_ms2_mz_met <- mat_raw_ms2_mz[i,]
    mat_raw_ms2_mz_met <- na.omit(mat_raw_ms2_mz_met)
    mat_raw_ms2_int_met <- mat_raw_ms2_int[i,]
    mat_raw_ms2_int_met <- na.omit(mat_raw_ms2_int_met)
    mat_raw_ms2_int_met <- as.numeric(mat_raw_ms2_int_met)
    mat_raw_ms2_int_met_hight <- which(mat_raw_ms2_int_met == max(mat_raw_ms2_int_met))
    mat_raw_ms2_10_hight <- as.numeric(mat_raw_ms2_int_met[mat_raw_ms2_int_met_hight])*(int_cut_level/100)
    mat_raw_ms2_mz_met_1 <- matrix(nrow = 1 ,ncol = length(mat_raw_ms2_int_met))
    ms2n <- 1
    ms2na <- 1
    for(ms2n in c(1:length(mat_raw_ms2_int_met))){
      if(mat_raw_ms2_int_met[ms2n] >= mat_raw_ms2_10_hight){
        mat_raw_ms2_mz_met_1[1,ms2na] <- mat_raw_ms2_mz_met[ms2n] 
        ms2na = ms2na + 1
      }else{
        
      }
    }
    mat_raw_ms2_mz_met <- na.omit(as.vector(mat_raw_ms2_mz_met_1))
    metabolite_id_result_table_1 <- matrix(0, nrow = length(fragment_list)+1 ,ncol = (reaction_num_1+1))
    
    a=1
    b=2
    c=3
    for(a in c(1:(reaction_num_1+1))){
      if(a==1){
        metabolite_id_result_table_1[1,a] <- paste0(mat_list_da[i,1])
      }else{
        if(reaction_num_1 == 1){
          metabolite_id_result_table_1[1,a] <- paste0(mat_list_da[i,1],"+",paste(na.omit(reaction_type_table_1[a-1]),collapse="+"))
        }else{
          metabolite_id_result_table_1[1,a] <- paste0(mat_list_da[i,1],"+",paste(na.omit(reaction_type_table_1[a-1,]),collapse="+"))
        }
      }
      for(b in c(1:length(fragment_list)+1)){
        for(c in c(1:length(mat_raw_ms2_mz_met))){
          ppm <- abs((round(as.numeric(mat_raw_ms2_mz_met[c]),4)-round(as.numeric(mat_pd_structure_ass[b-1,a]),4))/round(as.numeric(mat_raw_ms2_mz_met[c]),4))*1000000
          if(ppm <= ppm_value){
            metabolite_id_result_table_1[b,a] <- 1
          }else if(metabolite_id_result_table_1[b,a]==1){
            
          }else{
            metabolite_id_result_table_1[b,a] <- 0
          }
          
        }
      }
    }
    
    jj = 1
    jx = 1
    for(jj in c(1:(reaction_num_1))){
      if(reaction_num_1 == 1){
        reaction_posit <- which(meta_comb_list_f[,2] == paste(na.omit(reaction_type_table_1[reaction_num_1]),collapse="+"))
      }else{
        reaction_posit <- which(meta_comb_list_f[,2] == paste(na.omit(reaction_type_table_1[reaction_num_1,]),collapse="+"))
      }
      reaction_mz_vzlue <- round(as.numeric(meta_comb_list_f[reaction_posit,3]),4)-as.numeric(mat_raw_ms1_mz[1,1])
      
      for(jx in c(1:(length(fragment_list)))){
        if(as.numeric(metabolite_id_result_table_1[jx+1,1])==1 && as.numeric(metabolite_id_result_table_1[jx+1,jj+1])==1){
          
          if(reaction_mz_vzlue < 0){
            
            jxx = 1
            for(jxx in c(1:jx)){
              test_ppm <- abs((round(as.numeric(metabolite_id_result_table[jxx+1,1]),4) - (round(as.numeric(metabolite_id_result_table[jx+1,1]),4) + reaction_mz_vzlue))/(round(as.numeric(metabolite_id_result_table[jx+1,1]),4) + reaction_mz_vzlue))
              if(test_ppm <= ppm_value){
                metabolite_id_result_table_1[jx+1,jj+1] <- 0
              }else{
                next
              }
            }
            
          }else if(reaction_mz_vzlue > 0){
            
            jxy = jx
            for(jxy in c(jx:(length(fragment_list)))){
              if(round(as.numeric(metabolite_id_result_table[jxy+1,1]),4) == (round(as.numeric(metabolite_id_result_table[jx+1,1]),4) + reaction_mz_vzlue)){
                metabolite_id_result_table_1[jx+1,jj+1] <- 0
              }else{
                next
              }
            }
          }else{
            
          }
        }else{
          
        }
      }
    }
    metabolite_id_result_table <- cbind(metabolite_id_result_table,metabolite_id_result_table_1)
    
    metabolite_id_result_posit_1 <- matrix(nrow = 2 ,ncol = (reaction_num_1+1))
    metabolite_id_result_posit_1[1,] <- metabolite_id_result_table_1[1,]
    tt=1
    for(tt in c(1:(reaction_num_1+1))){
      if(tt == 1 && parent_site_uni == 1){
        frag_post <- which(metabolite_id_result_table_1[,tt] == 1)
        if(!(length(frag_post))){
          metabolite_id_result_posit_1[2,tt] <- 0
          next
        }else if(length(frag_post) == 1){
          nrow_frag_post = 1
        }else{
          nrow_frag_post <- nrow(mat_raw_ms2_pd_structure[frag_post-1,2:100])
        }
        frag_post_table <- mat_raw_ms2_pd_structure[frag_post-1,2:100]
        tts=1
        if(is.null(nrow_frag_post)){
          nrow_frag_post = 1
        }else{
          
        }
        
        if(nrow_frag_post>1){
          for(tts in c(1:(nrow_frag_post-1))){
            if(tts==1){
              overlap_frg_posit <- union(na.omit(frag_post_table[tts,]),na.omit(frag_post_table[tts+1,]))
            }else{
              overlap_frg_posit <- union(overlap_frg_posit,na.omit(frag_post_table[tts+1,]))
            }
          }
          
        }else{
          overlap_frg_posit <- na.omit(frag_post_table)
        }
        frag_post_v_crt <- paste(overlap_frg_posit,collapse="\t")
        metabolite_id_result_posit_1[2,tt] <- frag_post_v_crt
      }else{
        frag_post <- which(metabolite_id_result_table_1[,tt] == 1)
        if(!(length(frag_post))){
          metabolite_id_result_posit_1[2,tt] <- 0
          next
        }else if(length(frag_post) == 1){
          nrow_frag_post = 1
        }else{
          nrow_frag_post <- nrow(mat_raw_ms2_pd_structure[frag_post-1,2:100])
        }
        frag_post_table <- mat_raw_ms2_pd_structure[frag_post-1,2:100]
        tts=1
        if(is.null(nrow_frag_post)){
          nrow_frag_post = 1
        }else{
          
        }
        
        if(nrow_frag_post>1){
          for(tts in c(1:(nrow_frag_post-1))){
            if(tts==1){
              overlap_frg_posit <- intersect(na.omit(frag_post_table[tts,]),na.omit(frag_post_table[tts+1,]))
            }else{
              overlap_frg_posit <- intersect(overlap_frg_posit,na.omit(frag_post_table[tts+1,]))
            }
          }
          
        }else{
          overlap_frg_posit <- na.omit(frag_post_table)
        }
        frag_post_v_crt <- paste(overlap_frg_posit,collapse="\t")
        metabolite_id_result_posit_1[2,tt] <- frag_post_v_crt
      }
    }
    metabolite_id_result_posit <- cbind(metabolite_id_result_posit, metabolite_id_result_posit_1)
    metabolite_id_result_short_score <- matrix(nrow = reaction_num_1 + 1 ,ncol = 1)
    kksts <- reaction_num_1 + 1
    sst = 1
    for (sst in c(1:kksts)){
      if(length(unlist(strsplit(metabolite_id_result_posit_1[2,sst],"\t", fixed = TRUE)))==0 || as.numeric(unlist(strsplit(metabolite_id_result_posit_1[2,sst],"\t", fixed = TRUE)))==0){
        rec_score = 0
      }else{
        rec_score = 1
      }
      metabolite_id_result_short_score[sst,1] <- rec_score
    }
    metabolite_id_result_table_score[1,cc] <- as.character(mat_list_da[i,(cc+2)])
    metabolite_id_result_table_score_p[1,cc] <- as.character(mat_list_da[i,(cc+2)])
    metabolite_id_result_table_score[i,cc] <- sum(as.numeric(metabolite_id_result_short_score[,1]))
    metabolite_id_result_table_score_p[i,cc] <- sum(as.numeric(metabolite_id_result_short_score[,1]))/(reaction_num_1 + 1)
    
    #calcultion of cosine
    cosine_p <- vector("numeric", (reaction_num_1 + 1))
    cosine_data <- vector("numeric", (reaction_num_1 + 1))
    vi = 1
    for(vi in c(1:(reaction_num_1 + 1))){
      cosine_p[vi] <- 1
    }
    
    metabolite_id_result_posit_a1 <- t(metabolite_id_result_posit_1)
    vii = 1
    for(vii in c(1:(reaction_num_1 + 1))){
      sum_cosine_data <- sum(as.numeric(unlist(strsplit(as.character(metabolite_id_result_posit_a1[vii,2]),"\t", fixed = TRUE))))
      length_cosine_data <- length(as.numeric(unlist(strsplit(as.character(metabolite_id_result_posit_a1[vii,2]),"\t", fixed = TRUE))))
                                                          
       if(sum_cosine_data > 0 && length_cosine_data > 1){
          cosine_data[vii] <- 1
       }else{
          cosine_data[vii] <- 0
       }
    }
    metabolite_id_result_table_score_p_1[i,cc] <- round(as.numeric(cosine(cosine_p,cosine_data)),3)
  }

  assign(metabolite_id_result_posit_name, t(metabolite_id_result_posit))
  assign(metabolite_id_result_table_name, t(metabolite_id_result_table))
}

plot(pd_structure_smile ,atomcex=1, colbonds= as.numeric(unlist(strsplit(as.character(metabolite_id_result_posit_f_1[2,2]),"\t", fixed = TRUE))), print=FALSE)