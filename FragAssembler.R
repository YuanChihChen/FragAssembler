#==========================================================
# 1. install packages
#.libPaths("/usr/lib/R/site-library")
#options(repos = c(CRAN = "https://cran.r-project.org")) ##=>?���?????��??
#getOption("repos") ##=>?���?????��??

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager", dependencies = TRUE) ##=>修改??��?�dependencies
# BiocManager::install(version = "3.17",ask = FALSE)
# BiocManager::install("ChemmineR")

# install.packages("lsa")

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager",dependencies = TRUE) ##=>修改??��?�dependencies
# BiocManager::install("fmcsR") 


# BiocManager::install("ChemmineOB")
#==========================================================

#==========================================================
# 2. loading packages & function
library("ChemmineR")
library("lsa")
library("ChemmineOB")
library("fmcsR") 
options(show.error.locations = TRUE)
pascalTriangle_1 <- function(h) {
  lapply(h, function(i) choose(i, 1:i))
}
#==========================================================
args = commandArgs(trailingOnly=TRUE)
#==========================================================
# 3. loading tables
meta_pathy_list <-as.matrix(read.csv(paste0(args[7],"Biotransformation database.csv"), header=T, sep=",")) #input Biotransformation database table
meta_cand_list <-as.matrix(read.csv(paste0(args[7],"Biotransformation product candidates list_mehp.csv"), header=T, sep=","))  #input Biotransformation product candidates list table
msms_frag_data <- read.csv(paste0(args[7],"sfd.csv"), header=F, sep=",") #msms structure table
#==========================================================

#==========================================================
# 4. setting parameters
ppm_value <- strtoi(args[1])#5 
roop_num <- strtoi(args[2])#4
phase_i <- strtoi(args[3])#1
phase_ii <- strtoi(args[4])#0
mat_raw_charge <- args[5]#"CHARGE=1-"
int_cut_level <- strtoi(args[6])#0
#==========================================================


# 5. start fragassembler calculation
sss=1
for(sss in c(1:1)){
  ppm_data <- ppm_value

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
  if(phase_i == 0 && meta_pathy_list[i,3]==1){
    meta_pathy_list_new <- meta_pathy_list_new[-ss,]
    ss = ss-1
  }else if (phase_ii == 0 && meta_pathy_list[i,3]==2){
    meta_pathy_list_new <- meta_pathy_list_new[-ss,]
    ss = ss-1
  }
  ss = ss + 1
}
meta_pathy_list_num <- length(meta_pathy_list_new[,1])
meta_cand_list_num <- length(meta_cand_list[,1])
precursor_data <- meta_cand_list[1,]

#Generating the list of all possible metabolic reactions
j=1
i=1
full_list_num_all <- meta_pathy_list_num^roop_num
meta_comb_list <- matrix(nrow = full_list_num_all*2 ,ncol = 3)
meta_comb_list_count <- matrix(nrow = full_list_num_all*2 ,ncol = (roop_num+1))
colnames(meta_comb_list) = c("count", "pathway", "mz")
count_a <- 1
for(i in c(1:roop_num)){
  full_list_num = meta_pathy_list_num^i
  if(i==1){
    for(j in c(1:full_list_num)){
      meta_comb_list[count_a,1]<- count_a
      meta_comb_list[count_a,2]<- meta_pathy_list_new[j,1]
      meta_comb_list[count_a,3]<- round(as.numeric(precursor_data[2]),5)+round(as.numeric(meta_pathy_list_new[j,2]),5)
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

# rempve duplicated possible metabolic reactions
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
meta_comb_list_f <- matrix(nrow = meta_comb_list_num ,ncol = 3)
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

# Matching m/z candidates to all possible metabolic reactions
meta_cand_list_num <- length(meta_cand_list[,1])
meta_cand_id_result_mz_match <- matrix(nrow = meta_cand_list_num*100 ,ncol = 14)
meta_cand_id_result_pathway_length_mz <- matrix(nrow = meta_cand_list_num*100 ,ncol = meta_cand_list_num)
i=1
for(i in c(1:meta_cand_list_num)){
  meta_cand_id_result_pathway_length_mz[1,i] <- as.character(meta_cand_list[i,1])
}

mz_match_cc <- 1

i=1
j=1
for(i in c(1:meta_cand_list_num)){
  ddk <- 0
  ppm_old = "x"
  mz_match_cc_a <- 1
  for(j in c(1:meta_comb_list_f_num)){
    ppm <- abs(1000000*((as.numeric(meta_comb_list_f[j,3])-as.numeric(meta_cand_list[i,2]))/as.numeric(meta_comb_list_f[j,3])))
    ppm_1 <- 1000000*((as.numeric(meta_comb_list_f[j,3])-as.numeric(meta_cand_list[i,2]))/as.numeric(meta_comb_list_f[j,3]))
    if(ppm <= ppm_data && ppm_old < ppm){
      ddk = 10
      mz_match_cc = mz_match_cc - 1
      meta_cand_id_result_mz_match[mz_match_cc,1] <- meta_cand_list[i,1]
      meta_cand_id_result_mz_match[mz_match_cc,2] <- round(as.numeric(meta_comb_list_f[j,3]), digits=4)
      meta_cand_id_result_mz_match[mz_match_cc,3] <- "look right"
      meta_cand_id_result_mz_match[mz_match_cc,4:12] <- 0
      meta_cand_id_result_mz_match[mz_match_cc,13] <- meta_comb_list_f[j,2]
      meta_cand_id_result_mz_match[mz_match_cc,14] <- ppm_1
      meta_cand_id_result_pathway_length_mz[mz_match_cc_a+1,i] <- meta_comb_list_f[j,2]
      mz_match_cc_a = mz_match_cc_a + 1
      mz_match_cc = mz_match_cc + 1
      #ppm_old <- ppm
    }else if(ppm <= ppm_data && ppm_old == "x"){
      ddk = 10
      meta_cand_id_result_mz_match[mz_match_cc,1] <- meta_cand_list[i,1]
      meta_cand_id_result_mz_match[mz_match_cc,2] <- round(as.numeric(meta_comb_list_f[j,3]), digits=4)
      meta_cand_id_result_mz_match[mz_match_cc,3] <- "look right"
      meta_cand_id_result_mz_match[mz_match_cc,4:12] <- 0
      meta_cand_id_result_mz_match[mz_match_cc,13] <- meta_comb_list_f[j,2]
      meta_cand_id_result_mz_match[mz_match_cc,14] <- ppm_1
      meta_cand_id_result_pathway_length_mz[mz_match_cc_a+1,i] <- meta_comb_list_f[j,2]
      mz_match_cc_a = mz_match_cc_a + 1
      mz_match_cc = mz_match_cc + 1
      #ppm_old <- ppm
    }else{
    }
    
  }
  if(ddk == 0){

    meta_cand_id_result_mz_match[mz_match_cc,1:3] <- meta_cand_list[i,1:3]
    meta_cand_id_result_mz_match[mz_match_cc,13] <- "ND"
    meta_cand_id_result_mz_match[mz_match_cc,14] <- "ND"
    meta_cand_id_result_pathway_length_mz[mz_match_cc_a+1,i] <- "ND"
    mz_match_cc_a = mz_match_cc_a + 1
    mz_match_cc = mz_match_cc + 1
  }else{
    
  }
  
}
meta_cand_id_result_mz_match <- meta_cand_id_result_mz_match[complete.cases(meta_cand_id_result_mz_match[,1]),]
#generate the list for following ms/ms calculation
i=2
calculate_reaction_num_table <- matrix(nrow = (meta_cand_list_num-1) ,ncol = 1)
for(i in c(2:meta_cand_list_num)){
  calculate_reaction_num_table[(i-1),1] <- length(na.omit(meta_cand_id_result_pathway_length_mz[,i]))
}

calculate_reaction_num <- as.numeric(max(calculate_reaction_num_table[1:(meta_cand_list_num-1),1]))

precursor_num <- 1

mat_list_da <- matrix(nrow = meta_cand_list_num ,ncol = calculate_reaction_num+2)

mat_list_da[1,1]<- meta_cand_list[precursor_num,1]
mat_list_da[1,2]<- meta_cand_list[precursor_num,3]
i=1
for(i in c(1:calculate_reaction_num)){
  mat_list_da[1,i+2] <- "ND"
}

aa=1
input_code = 1 
for(aa in c(1:meta_cand_list_num)){
  if(aa == precursor_num){
    input_code = input_code + 1
  }else{ 
    mat_list_da[input_code,1]<- meta_cand_list[aa,1]
    mat_list_da[input_code,2]<- meta_cand_list[aa,3]

    j=1
    for(j in c(1:calculate_reaction_num)){
      if(is.na(meta_cand_id_result_pathway_length_mz[j+1,aa])){
        mat_list_da[input_code,j+2] <- "ND"
      }else{
        mat_list_da[input_code,j+2] <- meta_cand_id_result_pathway_length_mz[j+1,aa]
      } 
    }
    input_code = input_code + 1
  }
}



# assign structure of parent compound fragment
mat_list_da_num <- length(mat_list_da[,1])

msms_frag_data_num_org <- length(msms_frag_data[,1]) # msms data arrangement

msms_frag_data_1 <- matrix(nrow = msms_frag_data_num_org ,ncol = 22)

i=1
for(i in c(1:msms_frag_data_num_org)){
  n_count_split <- length(unlist(strsplit(as.character(msms_frag_data[i,1]),"\t")))
  msms_frag_data_1[i,1:n_count_split] <- unlist(strsplit(as.character(msms_frag_data[i,1]),"\t"))
}


isst = 1
for(isst in c(1:msms_frag_data_num_org)){
  test_line <- as.character(msms_frag_data_1[isst,1])
  test_line_1 <-strsplit(test_line,"VS2")
  test_line_2 <-strsplit(test_line,": ")
  if(test_line_1[[1]][1] == "Num Fragment "){
    msms_frag_data_cunt <- isst + 1
    break
  }else if(test_line_2[[1]][1] == "SMILES"){
    precursor_smiles <- test_line_2[[1]][2]
  }else{
    
  }
}
msms_frag_data_1 <- cbind(msms_frag_data_1[msms_frag_data_cunt:msms_frag_data_num_org,1],msms_frag_data_1[msms_frag_data_cunt:msms_frag_data_num_org,16])
msms_frag_data_num <- length(msms_frag_data_1[,1])

mat_raw_ms2_pd_structure <- matrix(nrow = msms_frag_data_num ,ncol = 100)

as(precursor_smiles, "SMIset")
pd_structure_smile <- smiles2sdf(precursor_smiles)

i=1
for(i in c(1:msms_frag_data_num)){
  as(msms_frag_data_1[i,2], "SMIset")
  frag_structure_smile <- smiles2sdf(msms_frag_data_1[i,2])
  test <- fmcs(pd_structure_smile, frag_structure_smile, au=0, bu=0)
  mcsa1 <- mcs1(test)
  frag_structure_code <- as.numeric(unlist(mcsa1[[2]][1]))
  frag_structure_length <- length(as.numeric(unlist(mcsa1[[2]][1])))
  mat_raw_ms2_pd_structure[i,1] <- as.numeric(msms_frag_data_1[i,1])
  mat_raw_ms2_pd_structure[i,2:(frag_structure_length+1)] <- frag_structure_code
}

fragment_mz_list <- matrix(nrow = msms_frag_data_num ,ncol = 1)
fragment_mz_list[,1] <- mat_raw_ms2_pd_structure[,1]

parent_mz <-round(as.numeric(meta_cand_list[1,2]),4)

# input experiment msms data
mat_raw_ms1_mz <- matrix(nrow = mat_list_da_num ,ncol = 1000)
mat_raw_ms1_int <- matrix(nrow = mat_list_da_num ,ncol = 1000)
mat_raw_ms2_mz <- matrix(nrow = mat_list_da_num ,ncol = 1000)
mat_raw_ms2_int <- matrix(nrow = mat_list_da_num ,ncol = 1000)
data_input_num_ms1 <- 1000
data_input_num_ms2 <- 1000

i=1
for(i in c(1:mat_list_da_num)){
  mat_pathway <- mat_list_da[i,2]
  mat_list_da_n <- mat_list_da[i,1]
  mat_raw <- readLines(mat_pathway) 
  mat_raw_num <- length(mat_raw)
  x=1
  vvv <- mat_raw_num + 1
  ooo = 0
  for(x in c(1:mat_raw_num)){
    xxx <- unlist(strsplit(as.character(mat_raw[(vvv-x)]),"\t"))
    if(is.na(xxx[2]) & ooo == 0){
      mat_raw <- mat_raw[-(vvv-x)]
      mat_raw_num <- length(mat_raw)
    }else if(!is.na(xxx[2]) & ooo == 0){
      ooo <- 1
    }
  }
  
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
    mat_list_damsc2 <- unlist(strsplit(as.character(mat_raw[j-1]),":"))
    
    if(j >= 3 && mat_list_damsc1[1] == "MSTYPE: MS1"){
      record_ms1 <- 1
      ms1_number_mat <- as.numeric(mat_list_damsc2[2])
    }else if(j >= 3 && mat_list_damsc[1] == "MSTYPE: MS2"){
      record_ms1 <- 0
    }else if(j >= 3 && mat_list_damsc1[1] == "MSTYPE: MS2"){
      ms2_number_mat <- as.numeric(mat_list_damsc2[2])
      record_ms2 <- 1
    }else{
      
    }
    
    if(record_ms1 == 1){
      mat_raw_ms1_rawvalue <- unlist(strsplit(as.character(mat_raw[j]),"\t"))
      if(ms1_number_mat > data_input_num_ms1 & cord_ms1 == 1){
        mat_raw_ms1_mz_1 <- matrix(nrow = mat_list_da_num ,ncol = ms1_number_mat)
        mat_raw_ms1_int_1 <- matrix(nrow = mat_list_da_num ,ncol = ms1_number_mat)
        mat_raw_ms1_mz <- cbind(mat_raw_ms1_mz,mat_raw_ms1_mz_1)
        mat_raw_ms1_int <- cbind(mat_raw_ms1_int,mat_raw_ms1_int_1)
        data_input_num_ms1 <- data_input_num_ms1 + ms1_number_mat
        mat_raw_ms1_mz[i,cord_ms1] <- mat_raw_ms1_rawvalue[1] 
        mat_raw_ms1_int[i,cord_ms1] <- mat_raw_ms1_rawvalue[2]
      }else{
        mat_raw_ms1_mz[i,cord_ms1] <- mat_raw_ms1_rawvalue[1] 
        mat_raw_ms1_int[i,cord_ms1] <- mat_raw_ms1_rawvalue[2]
      }
      cord_ms1 = cord_ms1 +1
    }else if(record_ms2 == 1){
      mat_raw_ms2_rawvalue <- unlist(strsplit(as.character(mat_raw[j]),"\t"))
      if(ms2_number_mat > data_input_num_ms2 & cord_ms2 == 1){
        mat_raw_ms2_mz_1 <- matrix(nrow = mat_list_da_num ,ncol = ms2_number_mat)
        mat_raw_ms2_int_1 <- matrix(nrow = mat_list_da_num ,ncol = ms2_number_mat)
        mat_raw_ms2_mz <- cbind(mat_raw_ms2_mz,mat_raw_ms2_mz_1)
        mat_raw_ms2_int <- cbind(mat_raw_ms2_int,mat_raw_ms2_int_1)
        data_input_num_ms2 <- data_input_num_ms2 + ms2_number_mat       
        mat_raw_ms2_mz[i,cord_ms2] <- mat_raw_ms2_rawvalue[1] 
        mat_raw_ms2_int[i,cord_ms2] <- mat_raw_ms2_rawvalue[2]
      }else{
        mat_raw_ms2_mz[i,cord_ms2] <- mat_raw_ms2_rawvalue[1] 
        mat_raw_ms2_int[i,cord_ms2] <- mat_raw_ms2_rawvalue[2]
      }
      cord_ms2 = cord_ms2 +1
    }else{
      
    }
  }
}
pd_ms2_num <- length(na.omit(mat_raw_ms2_mz[1,]))

msms_data_pre_exp <- mat_raw_ms2_mz[1,]
msms_data_pre_exp <- na.omit(msms_data_pre_exp)
msms_data_pre_exp_num <- length(msms_data_pre_exp)
msms_data_pre_exp_match <- matrix(nrow = msms_frag_data_num ,ncol = 3)
msms_data_pre_exp_match[,1] <- mat_raw_ms2_pd_structure[,1]

ppm_mem <- 10

i=1
j=1
for(i in c(1:msms_frag_data_num)){
  
  for(j in c(1:msms_data_pre_exp_num)){
    ppm_check <- abs((round(as.numeric(msms_data_pre_exp[j]),4)-round(as.numeric(mat_raw_ms2_pd_structure[i,1]),4))/round(as.numeric(mat_raw_ms2_pd_structure[i,1]),4))*1000000
    if(ppm_check<=3 && ppm_mem == 10){
      msms_data_pre_exp_match[i,2] <- j
      msms_data_pre_exp_match[i,3] <- ppm_check
      ppm_mem <- ppm_check
    }else if(ppm_check<=3 && ppm_mem != 10 && ppm_mem >= ppm_check){
      msms_data_pre_exp_match[i,2] <- j
      msms_data_pre_exp_match[i,3] <- ppm_check
      ppm_mem <- ppm_check
    }else{
      
    }  
  }
}

fragment_list <- as.numeric(t(msms_data_pre_exp_match[,2]))

# idntification of fragment signatures
metabolite_id_result_table_score <- matrix(nrow = mat_list_da_num ,ncol = calculate_reaction_num)
metabolite_id_result_table_score_p <- matrix(nrow = mat_list_da_num ,ncol = calculate_reaction_num)

metabolite_id_result_table_score_p_1 <- matrix(nrow = mat_list_da_num ,ncol = calculate_reaction_num)

cc=1
for(cc in c(1:calculate_reaction_num)){
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
    if(as.character(mat_list_da[i,(cc+2)]) == "ND"){
      next
    }else{
      
    }
 
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
      reaction_mz_vzlue <- round(as.numeric(meta_comb_list_f[reaction_posit,3]),4)-parent_mz
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
    metabolite_id_result_posit <- cbind(metabolite_id_result_posit, metabolite_id_result_posit_1)
    metabolite_id_result_short_score <- matrix(nrow = reaction_num_1 + 1 ,ncol = 1)
    kksts <- reaction_num_1 + 1
    sst = 1
    for (sst in c(1:kksts)){
      if((length(unlist(strsplit(metabolite_id_result_posit_1[2,sst],"\t", fixed = TRUE)))==0) || any((as.numeric(unlist(strsplit(metabolite_id_result_posit_1[2,sst],"\t", fixed = TRUE)))==0))){   # <- ??�自行修?��(??��?�any，�?�????�詢???)
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
    vi = 2
    for(vi in c(2:(reaction_num_1 + 1))){
      cosine_p[vi] <- 1
    }
    
    metabolite_id_result_posit_a1 <- t(metabolite_id_result_posit_1)
    vii = 2
    for(vii in c(2:(reaction_num_1 + 1))){
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
if(nrow(metabolite_id_result_table_score_p_1) == 2){
  metabolite_id_result_table_score_p_2 <- t(as.matrix(metabolite_id_result_table_score_p_1[-1,]))
}else{
  metabolite_id_result_table_score_p_2 <- metabolite_id_result_table_score_p_1[-1,]
}

meta_cand_id_result_mz_match <- meta_cand_id_result_mz_match[,-3:-12]
rownames(metabolite_id_result_table_score_p_2) <- mat_list_da[2:mat_list_da_num,1]
colnames(metabolite_id_result_table_score_p_2) <- c(1:calculate_reaction_num)
colnames(meta_cand_id_result_mz_match) = c("name", "m/z",	"metabolic reactions",	"mass error")

output_results_table_total <- matrix(nrow = (mat_list_da_num-1) ,ncol = 3)
if(nrow(mat_list_da) == 2){
  output_results_table_total_reaction <- t(as.matrix(mat_list_da[2:mat_list_da_num, 3:(calculate_reaction_num+2)]))
}else{
  output_results_table_total_reaction <- mat_list_da[2:mat_list_da_num, 3:(calculate_reaction_num+2)]
}

colnames(output_results_table_total_reaction) <- paste("reaction_", 1:calculate_reaction_num, sep = "")

output_results_table_total[,1] <- paste("M_", 1:(mat_list_da_num-1), sep = "")
output_results_table_total[,2] <- meta_cand_list[2:mat_list_da_num,2]
output_results_table_total[,3] <- mat_list_da[2:mat_list_da_num,2]
colnames(output_results_table_total) <- c("Name", "m/z", "MS/MS file name")
output_results_table_total <- cbind(output_results_table_total,output_results_table_total_reaction)

# generate result table

result_list <- list()
i=1
for(i in c(1:(mat_list_da_num-1))){
  mat_result <- matrix(nrow = 1 ,ncol = 3)
  result_list[[i]] <- mat_result
  result_list[[i]][1,1] <- "reaction number"
  result_list[[i]][1,2] <- "name"
  result_list[[i]][1,3] <- "site"
}

i=1
for(i in c(1:calculate_reaction_num)){
  result_transfer <- get(paste("metabolite_id_result_posit_f",i, sep = "_"))
  result_transfer_num <- length(result_transfer[,1])
  if(length(result_transfer)==2){
    break
  }else{
    
  }
  
  j=2
  for(j in c(2:result_transfer_num)){
    result_transfer_cand_first <- unlist(strsplit(as.character(result_transfer[j,1]),"+", fixed = TRUE))
    if(length(result_transfer_cand_first)==1){
      next
    }else{
      metabolite_count <- (which(meta_cand_list[,1] == result_transfer_cand_first[1]))-1
      metabolite_count_name <- paste0("M_", metabolite_count)
      k=1
      result_transfer_cand_first_name <- "name"
      for(k in c(1:length(result_transfer_cand_first))){
        if(k==1){
          result_transfer_cand_first_name <- metabolite_count_name
        }else{
          result_transfer_cand_first_name <- paste0(result_transfer_cand_first_name, "+", result_transfer_cand_first[k])
        }
        
      }
    }
    result_transfer_cand <- unlist(strsplit(as.character(result_transfer[j,]),"+", fixed = TRUE))
    result_transfer_posit <- which(meta_cand_list[,1] == result_transfer_cand[1])
    result_transfer_cand_1 <- matrix(nrow = 1 ,ncol = 3)
    result_transfer_cand_1 <- c(i, result_transfer_cand_first_name, result_transfer[j,2])
    result_transfer_cand_2 <- rbind(result_list[[(result_transfer_posit-1)]],result_transfer_cand_1)
    rownames(result_transfer_cand_2) <- NULL
    result_list[[(result_transfer_posit-1)]]<-result_transfer_cand_2
  }
}
compare_name_table <- matrix(nrow = (mat_list_da_num-1) ,ncol = 2)
compare_name_table[,2] <- output_results_table_total[,1]

}

# 6. output tables
write.csv(output_results_table_total,paste0(args[7],"FragAssembler_metabolites_information.csv") , row.names = FALSE) #export result, in "Documents"  "Assigning_biotransformation_reaction_mass.csv"
write.csv(metabolite_id_result_table_score_p_2,paste0(args[7],"FragAssembler_identification_score.csv") , row.names = FALSE) #export result, in "Documents"  "metabolite_id_result_table_score_p_1.csv"
write.SDF(pd_structure_smile,paste0(args[7],"pd_structure_smile.sdf"))
# 7. code for draw the FragAssembler results
# plot(pd_structure_smile ,atomcex=1, colbonds= as.numeric(unlist(strsplit(as.character(metabolite_id_result_posit_f_1[2,2]),"\t", fixed = TRUE))), print=FALSE)
output_file_name_table <- matrix(nrow = (mat_list_da_num-1) ,ncol = 2)
i=1
for(i in c(1:(mat_list_da_num-1))){
  outpute_file <- result_list[[i]]
  colnames(outpute_file) = c("reaction number","name", "position")
  outpute_file <- outpute_file[-1,]
  metabolite_id_result_posit_name_1 <-  paste0("metabolite_id_result_posit_f_",i)
  write.csv(outpute_file,paste0(args[7],paste0("/plot/",metabolite_id_result_posit_name_1,".csv")), row.names = FALSE)
}


output_results_table_list_each <- list()
i=1
for(i in c(1:(mat_list_da_num-1))){
  output_results_table_for_plot <- matrix(nrow = length(na.omit(meta_cand_id_result_pathway_length_mz[,(i+1)]))-1 ,ncol = 6)
  colnames(output_results_table_for_plot) = c("Reaction number","Metabolic reaction","Theoretical m/z","Measured m/z","Mass error (ppm)","Confidence_score")
  output_results_table_list_each[[i]] <- output_results_table_for_plot
  output_results_table_list_each[[i]][,1] <- c(1:(length(na.omit(meta_cand_id_result_pathway_length_mz[,(i+1)]))-1))
  output_results_table_list_each[[i]][,2] <- meta_cand_id_result_pathway_length_mz[2:(length(na.omit(meta_cand_id_result_pathway_length_mz[,(i+1)]))),(i+1)]
  output_results_table_list_each[[i]][,4] <- round(as.numeric(meta_cand_list[(i+1),2]),digits=4)
  j=1
  for(j in c(1:length(meta_cand_id_result_mz_match[,1]))){
    if(mat_list_da[(i+1),1]==meta_cand_id_result_mz_match[j,1]){
      output_results_table_list_each[[i]][,3] <- round(as.numeric(meta_cand_id_result_mz_match[j,2]),digits=4)
      output_results_table_list_each[[i]][,5] <- round(as.numeric(meta_cand_id_result_mz_match[j,4]),digits=2)
    }else{
      
    }
  }
  output_results_table_list_each[[i]][,6] <- metabolite_id_result_table_score_p_2[i,1:(length(na.omit(meta_cand_id_result_pathway_length_mz[,(i+1)]))-1)]
  write.csv(output_results_table_list_each[[i]],paste0(args[7],paste0("/table/metabolite_id_result_posit_show_",i,".csv")), row.names = FALSE)
  outpute_file_name <- unlist(strsplit(as.character(meta_cand_list[i+1,1]),".mat", fixed = TRUE))
  outpute_file_name <- chartr(".", "_", outpute_file_name)
  output_file_name_table[i,1] <- paste0(outpute_file_name[1],".csv")
  output_file_name_table[i,2] <- paste0("metabolite_id_result_posit_show_",i,".csv")
  compare_name_table[i,1]  <- paste0("metabolite_id_result_posit_show_",i,".csv")
  colnames(output_file_name_table) <- c("new","old")
}
write.csv(output_file_name_table,paste0(args[7],"output_file_name_table.csv"), row.names = FALSE)

write.csv(compare_name_table,paste0(args[7],"compare_name_table.csv"))
