meta_pathy_list <-as.matrix(read.csv(file.choose(), header=T, sep=",")) #input Biotransformation database table
meta_cand_list <-as.matrix(read.csv(file.choose(), header=T, sep=","))  #input Biotransformation product candidates list table
roop_num <- 4
phase <- 1
ppm_data <- 10

check.integer <- function(N){
  !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
}
check.integer_no <- function(N){
  grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
}
meta_pathy_list_num <- length(meta_pathy_list[,1])
i=1
kk <- 1
meta_pathy_list_new <- meta_pathy_list
ss=1
for(i in c(1:meta_pathy_list_num)){
  if(as.numeric(meta_pathy_list[i,3])>phase){
    meta_pathy_list_new <- meta_pathy_list_new[-ss,]
  }else{
    ss=ss+1
  }
}
meta_pathy_list_num <- length(meta_pathy_list_new[,1])
meta_cand_list_num <- length(meta_cand_list[,1])
precursor_data <- meta_cand_list[1,]

#generate all possible biotransformation
j=1
i=1
full_list_num_all <- meta_pathy_list_num^roop_num
meta_comb_list <- matrix(nrow = full_list_num_all*2 ,ncol = 12)
meta_comb_list_count <- matrix(nrow = full_list_num_all*2 ,ncol = (roop_num+1))
colnames(meta_comb_list) = c("count", "pathway", "mz", "C", "H", "N", "O", "P", "S", "F", "Br", "Cl")
count_a <- 1
for(i in c(1:roop_num)){
  full_list_num = meta_pathy_list_num^i
  if(i==1){
    for(j in c(1:full_list_num)){
      meta_comb_list[count_a,1]<- count_a
      meta_comb_list[count_a,2]<- meta_pathy_list_new[j,1]
      meta_comb_list[count_a,3]<- round(as.numeric(precursor_data[2]),5)+round(as.numeric(meta_pathy_list_new[j,2]),5)
      meta_comb_list[count_a,4]<- round(as.numeric(precursor_data[4]),0)+round(as.numeric(meta_pathy_list_new[j,4]),0)
      meta_comb_list[count_a,5]<- round(as.numeric(precursor_data[5]),0)+round(as.numeric(meta_pathy_list_new[j,5]),0)
      meta_comb_list[count_a,6]<- round(as.numeric(precursor_data[6]),0)+round(as.numeric(meta_pathy_list_new[j,6]),0)
      meta_comb_list[count_a,7]<- round(as.numeric(precursor_data[7]),0)+round(as.numeric(meta_pathy_list_new[j,7]),0)
      meta_comb_list[count_a,8]<- round(as.numeric(precursor_data[8]),0)+round(as.numeric(meta_pathy_list_new[j,8]),0)
      meta_comb_list[count_a,9]<- round(as.numeric(precursor_data[9]),0)+round(as.numeric(meta_pathy_list_new[j,9]),0)
      meta_comb_list[count_a,10]<- round(as.numeric(precursor_data[10]),0)+round(as.numeric(meta_pathy_list_new[j,10]),0)
      meta_comb_list[count_a,11]<- round(as.numeric(precursor_data[11]),0)+round(as.numeric(meta_pathy_list_new[j,11]),0)
      meta_comb_list[count_a,12]<- round(as.numeric(precursor_data[12]),0)+round(as.numeric(meta_pathy_list_new[j,12]),0)
      meta_comb_list_count[count_a,1]<- count_a
      meta_comb_list_count[count_a,2]<- j
      count_a = count_a + 1
    }
    fix_num_b2 <- meta_pathy_list_num
    full_list_num_old <- meta_pathy_list_num
    
  }else{
    fix_num_a1 <- 1
    fix_num_a2 <- meta_pathy_list_num
    fix_num_b1 <- fix_num_b2 + 1
    fix_num_b2_old <- fix_num_b2
    
    for(j in c(1:full_list_num)){
      
      if(fix_num_a1 == (fix_num_a2+1)){
        fix_num_a1 <- 1
        fix_num_b1 = fix_num_b1 + 1
      }else{
        
      }
      meta_comb_list[count_a,1]<- count_a
      meta_comb_list[count_a,2]<- paste0(meta_comb_list[(fix_num_b1-full_list_num_old),2],"+",meta_pathy_list_new[fix_num_a1,1])
      meta_comb_list[count_a,3]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),3]),5)+round(as.numeric(meta_pathy_list_new[fix_num_a1,2]),5)
      meta_comb_list[count_a,4]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),4]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,4]),0)
      meta_comb_list[count_a,5]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),5]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,5]),0)
      meta_comb_list[count_a,6]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),6]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,6]),0)
      meta_comb_list[count_a,7]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),7]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,7]),0)
      meta_comb_list[count_a,8]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),8]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,8]),0)
      meta_comb_list[count_a,9]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),9]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,9]),0)
      meta_comb_list[count_a,10]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),10]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,10]),0)
      meta_comb_list[count_a,11]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),11]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,11]),0)
      meta_comb_list[count_a,12]<- round(as.numeric(meta_comb_list[(fix_num_b1-full_list_num_old),12]),0)+round(as.numeric(meta_pathy_list_new[fix_num_a1,12]),0)
      meta_comb_list_count[count_a,1]<- count_a
      if(i==2){
        meta_comb_list_count[count_a,2] <- meta_comb_list_count[(fix_num_b1-full_list_num_old),2]
      }else{
        meta_comb_list_count[count_a,2:i] <- meta_comb_list_count[(fix_num_b1-full_list_num_old),2:i]
      }
      meta_comb_list_count[count_a,i+1] <- fix_num_a1
      count_a = count_a + 1
      fix_num_a1 = fix_num_a1 + 1
    }
    fix_num_b2 <- fix_num_b2_old + full_list_num
    full_list_num_old <- full_list_num
  }
}

meta_comb_list <- meta_comb_list[complete.cases(meta_comb_list[,1]),]
meta_comb_list_count <- meta_comb_list_count[complete.cases(meta_comb_list_count[,1]),]
meta_comb_list_num <- length(meta_comb_list[,1])

#rempve duplicated
meta_comb_list_count[is.na(meta_comb_list_count)] <- 0
meta_comb_list_count_1 <- matrix(nrow = roop_num ,ncol = meta_comb_list_num)
colnames(meta_comb_list_count_1) = meta_comb_list_count[,1]
meta_comb_list_count_1 <- t(meta_comb_list_count[,2:(roop_num+1)])
meta_comb_list_count_2 <- matrix(nrow = roop_num ,ncol = meta_comb_list_num)
i=1
for(i in c(1:meta_comb_list_num)){
  meta_comb_list_count_2[,i] <- meta_comb_list_count_1[order(meta_comb_list_count_1[,i],decreasing=T),i]
}
meta_comb_list_count_2 <- t(meta_comb_list_count_2)
rownames(meta_comb_list_count_2)<- meta_comb_list_count[,1]
meta_comb_list_count_3 <- meta_comb_list_count_2[!duplicated(meta_comb_list_count_2),]
meta_comb_list_count_end <- row.names(meta_comb_list_count_3)
meta_comb_list_f <- matrix(nrow = meta_comb_list_num ,ncol = 12)
i=1
kksy <-1
for(i in c(1:meta_comb_list_num)){
  if(meta_comb_list[i,1]==meta_comb_list_count_end[kksy]){
    meta_comb_list_f[kksy,] <- meta_comb_list[i,]
    kksy = kksy + 1
  }else{
    
  }
}
meta_comb_list_f <- meta_comb_list_f[complete.cases(meta_comb_list_f[,1]),]
meta_comb_list_f_num <- length(meta_comb_list_f[,1])

#match biotransformation database to biotransformation product candidates
meta_cand_list_num <- length(meta_cand_list[,1])
meta_cand_id_result_mole_match <- matrix(nrow = meta_cand_list_num*100 ,ncol = 14)
meta_cand_id_result_mz_match <- matrix(nrow = meta_cand_list_num*100 ,ncol = 14)
mole_match_cc <- 1
mz_match_cc <- 1

i=1
j=1
for(i in c(1:meta_cand_list_num)){
  dou_num <- (round(as.numeric(meta_cand_list[i,4]),digits=0)+1)+((round(as.numeric(meta_cand_list[i,6]),digits=0)+round(as.numeric(meta_cand_list[i,8]),digits=0))/2)-((round(as.numeric(meta_cand_list[i,5]),digits=0)+round(as.numeric(meta_cand_list[i,10]),digits=0)+round(as.numeric(meta_cand_list[i,11]),digits=0)+round(as.numeric(meta_cand_list[i,12]),digits=0))/2)
  ddc <- 0
  ddk <- 0
    ppm_old = "x"
    for(j in c(1:meta_comb_list_f_num)){
      #if(meta_comb_list_f[j,4]>=0 && meta_comb_list_f[j,5]>=0 && meta_comb_list_f[j,6]>=0 && meta_comb_list_f[j,7]>=0 && meta_comb_list_f[j,8]>=0 && meta_comb_list_f[j,9]>=0 && meta_comb_list_f[j,10]>=0 && meta_comb_list_f[j,11]>=0 && meta_comb_list_f[j,12]>=0){
      
      #}else{
        #next
      #}
      if(meta_comb_list_f[j,4]==meta_cand_list[i,4] && meta_comb_list_f[j,5]==meta_cand_list[i,5] && meta_comb_list_f[j,6]==meta_cand_list[i,6] && meta_comb_list_f[j,7]==meta_cand_list[i,7] && meta_comb_list_f[j,8]==meta_cand_list[i,8] && meta_comb_list_f[j,9]==meta_cand_list[i,9] && meta_comb_list_f[j,10]==meta_cand_list[i,10] && meta_comb_list_f[j,11]==meta_cand_list[i,11] && meta_comb_list_f[j,12]==meta_cand_list[i,12] && check.integer(dou_num)){
        ddc = 1
        ppm <- abs(100000*((round(as.numeric(meta_comb_list_f[j,3]), digits=4)-round(as.numeric(meta_cand_list[i,2]), digits=4))/round(as.numeric(meta_cand_list[i,2]), digits=4)))
        meta_cand_id_result_mole_match[mole_match_cc,1:12] <- meta_cand_list[i,1:12]
        meta_cand_id_result_mole_match[mole_match_cc,13] <- meta_comb_list_f[j,2]
        meta_cand_id_result_mole_match[mole_match_cc,14] <- ppm
        ppm <- "x"
        mole_match_cc = mole_match_cc + 1
      }else if(check.integer_no(dou_num) && j == meta_comb_list_f_num){
        ddc = 2
        meta_cand_id_result_mole_match[mole_match_cc,1:12] <- meta_cand_list[i,1:12]
        meta_cand_id_result_mole_match[mole_match_cc,13] <- "in source fragmentation"
        meta_cand_id_result_mole_match[mole_match_cc,14] <- "in source fragmentation"
        mole_match_cc = mole_match_cc + 1
      }else{
      }
        ppm <- abs(1000000*((round(as.numeric(meta_comb_list_f[j,3]), digits=4)-round(as.numeric(meta_cand_list[i,2]), digits=4))/round(as.numeric(meta_cand_list[i,2]), digits=4)))
        if(ppm <= ppm_data && ppm_old < ppm){
          ddk = 10
          mz_match_cc = mz_match_cc - 1
          meta_cand_id_result_mz_match[mz_match_cc,1:2] <- meta_cand_list[i,1:2]
          meta_cand_id_result_mz_match[mz_match_cc,3] <- "look right"
          meta_cand_id_result_mz_match[mz_match_cc,4:12] <- meta_comb_list_f[j,4:12]
          meta_cand_id_result_mz_match[mz_match_cc,13] <- meta_comb_list_f[j,2]
          meta_cand_id_result_mz_match[mz_match_cc,14] <- ppm
          mz_match_cc = mz_match_cc + 1
          #ppm_old <- ppm
        }else if(ppm <= ppm_data && ppm_old == "x"){
          ddk = 10
          meta_cand_id_result_mz_match[mz_match_cc,1:2] <- meta_cand_list[i,1:2]
          meta_cand_id_result_mz_match[mz_match_cc,3] <- "look right"
          meta_cand_id_result_mz_match[mz_match_cc,4:12] <- meta_comb_list_f[j,4:12]
          meta_cand_id_result_mz_match[mz_match_cc,13] <- meta_comb_list_f[j,2]
          meta_cand_id_result_mz_match[mz_match_cc,14] <- ppm
          mz_match_cc = mz_match_cc + 1
          #ppm_old <- ppm
        }else{
        }

    }
    if(ddc == 0 && ddk != 0){
      meta_cand_id_result_mole_match[mole_match_cc,1:12] <- meta_cand_list[i,1:12]
      meta_cand_id_result_mole_match[mole_match_cc,13] <- "ND"
      meta_cand_id_result_mole_match[mole_match_cc,14] <- "ND"
      mole_match_cc = mole_match_cc + 1
      
    }else if(ddk == 0 && ddc != 0){
      meta_cand_id_result_mz_match[mz_match_cc,1:12] <- meta_cand_list[i,1:12]
      meta_cand_id_result_mz_match[mz_match_cc,13] <- "ND"
      meta_cand_id_result_mz_match[mz_match_cc,14] <- "ND"
      mz_match_cc = mz_match_cc + 1
    }else if(ddk == 0 && ddc == 0){
      meta_cand_id_result_mole_match[mole_match_cc,1:12] <- meta_cand_list[i,1:12]
      meta_cand_id_result_mole_match[mole_match_cc,13] <- "ND"
      meta_cand_id_result_mole_match[mole_match_cc,14] <- "ND"
      mole_match_cc = mole_match_cc + 1
      meta_cand_id_result_mz_match[mz_match_cc,1:12] <- meta_cand_list[i,1:12]
      meta_cand_id_result_mz_match[mz_match_cc,13] <- "ND"
      meta_cand_id_result_mz_match[mz_match_cc,14] <- "ND"
      mz_match_cc = mz_match_cc + 1
    }
    
}
meta_cand_id_result_mole_match <- meta_cand_id_result_mole_match[complete.cases(meta_cand_id_result_mole_match[,1]),]
meta_cand_id_result_mz_match <- meta_cand_id_result_mz_match[complete.cases(meta_cand_id_result_mz_match[,1]),]
colnames(meta_cand_id_result_mole_match) = c("name", "m/z",	"molecular formula",	"C",	"H",	"N",	"O",	"P",	"S",	"F",	"Br",	"Cl",	"biotransformation pathway",	"mass error")
colnames(meta_cand_id_result_mz_match) = c("name", "m/z",	"molecular formula",	"C",	"H",	"N",	"O",	"P",	"S",	"F",	"Br",	"Cl",	"biotransformation pathway",	"mass error")
write.csv(meta_cand_id_result_mole_match, "Assigning_biotransformation_reaction_formula.csv", row.names = FALSE)
write.csv(meta_cand_id_result_mz_match, "Assigning_biotransformation_reaction_mass.csv", row.names = FALSE) #export result, in "Documents"
