#source("https://bioconductor.org/biocLite.R")
#biocLite("kebabs")
library(Biostrings)
library(stringr)
library(kebabs)
library(ShortRead)
library(Matrix)
library(parallel)
library(lme4)
#library(data.table)
#library(stringdist)
library(stringi)
library(parallel)
#library(pbapply)


setwd("C:/Users/liutongtong/Desktop/project1/1000x10/generate")
set.seed(6)
chr4 <- readDNAStringSet('na12878.chr4.fasta')
length_of_chr4 <- chr4@ranges@width
name_chr4 <- chr4@ranges@NAMES


#get 20M sequence from na12878 chr4
#length_of_subseq <- 2000
print("The first step is getting 20M subsequence of na12878")
length_of_subseq <- ceiling(20/181*chr4@ranges@width)
startposi_of_subseq <- 666
endposi_of_subseq <- startposi_of_subseq+length_of_subseq-1
one_subseq <- DNAStringSet(chr4[[1]][startposi_of_subseq:endposi_of_subseq])
one_subseq@ranges@NAMES <- name_chr4
print("The first step has finished.")
print("####################################################################")


#the probability of mutation equals to 0.001, get the other sequence
print("The second step is generating SNP")
print("The rate of mutation is 1/650")
p_of_mutation <- 1/650
num_of_mutation <- ceiling(one_subseq@ranges@width*p_of_mutation)
loci_of_mutation <- matrix(sample(c(TRUE,FALSE),replace=T,size = one_subseq@ranges@width,prob = c(1/650,649/650)),nrow = length(one_subseq),ncol=width(one_subseq), byrow=TRUE)
loci <- which(loci_of_mutation[1,] %in% TRUE)
orig_base <- DNAStringSet(one_subseq[[1]][loci])
mutation_base <- genRandBioSeqs("DNA", 1, rowSums(loci_of_mutation), biostring=TRUE,seed = 6)
print("It Is Checking Bases. And it may be takes long time. Hope you are patient.")
A <- DNAStringSet("A")
#for (i in 1:mutation_base@ranges@width) {
#  if(mutation_base[[1]][i]==orig_base[[1]][i])
#    if(mutation_base[[1]][i]==A)
#      mutation_base[[1]][i] <- "T"
#    else
#      mutation_base[[1]][i] <- "A"
#}
check_base <- function(var1){
  if(mutation_base[[1]][var1]==orig_base[[1]][var1]){
    if(mutation_base[[1]][var1]==A){
      var1 <- "T"
    }
    else
      var1 <- "A"
  }else{
    var1 <- as.character(mutation_base[[1]][var1])
  }
  return(var1)
}
v <- c(1:mutation_base@ranges@width)
start_time <- proc.time()
n_cores <- detectCores() - 1
#cl <- makeCluster(n_cores, type = "FORK")
cl <- makeCluster(n_cores,type = "PSOCK")
clusterEvalQ(cl,library(Biostrings))
clusterExport(cl, varlist=c("mutation_base","orig_base","A"))
mutation_base <- parLapply(cl,v,check_base)
stopCluster(cl)
end_time <- proc.time()
print("The time that Checking Bases spends is:")
print(end_time-start_time)
#start_time <- proc.time()
#mutation_base <- lapply(v,check_base)
#end_time <- proc.time()
#print(end_time-start_time)
mutation_base <- paste(unlist(mutation_base),collapse = "")
mutation_base <- DNAStringSet(mutation_base)
names(mutation_base) <- "mutation_base"
#orig_base <- DNAStringSet(orig_base[[1]])
#names(orig_base) <- "orig_base"
print("Checking Bases has finished!!!!!!!!!!!")
#n <- 0
#for (i in 1:mutation_base@ranges@width) {
#  if(mutation_base[[1]][i]==orig_base[[1]][i])
#    n <- n+1
#}
#print(n)
loci <- loci + startposi_of_subseq - 1
snp_info <- data.frame(loci,orig_base[[1]],mutation_base[[1]])
print("The snp informarion is outputting to a csv file.")
write.csv(snp_info,file = "snp_info.csv")
rm(snp_info)
print("The snp_info has been removed.")
two_subseq <- replaceLetterAt(one_subseq,loci_of_mutation,mutation_base)
#writeXStringSet(one_subseq, filepath ="na20M.fasta", append=FALSE,compress=FALSE, compression_level=NA, format="fasta")
print("The second step has finished.")
print("####################################################################")


print("The third step is generating the fragment matrix, and it will take much time.")
print("The matrix is 1000x10, that means 1000 times experiments and getting 10 20K-segments once.")
num_of_reads <- 10
#num_of_lab <- 100
num_of_lab <- 1000
#length_of_reads <- 10
length_of_reads <- 20*1024
en1 <- length_of_subseq-length_of_reads+1
#one_subseq <- genRandBioSeqs("DNA", 1, length_of_subseq, biostring=TRUE,seed = 6)
#two_subseq <- genRandBioSeqs("DNA", 1, length_of_subseq, biostring=TRUE,seed = 8)
reads_matrix <- replicate(num_of_lab,{a <- DNAStringSet()})

generate_reads <- function(){
  sel_allele <- sample(1:2,size = 1,replace = T)
  loci_r <- sample(1:en1,size = 1,replace = T)
  en2 <- loci_r-1+length_of_reads
  if(sel_allele==1) {
    tmp <- DNAStringSet(one_subseq[[1]][loci_r:en2])
  } else {
    tmp <- DNAStringSet(two_subseq[[1]][loci_r:en2])
  }
  tmp@ranges@NAMES <- name_chr4
  return(tmp)
}

generate_reads_matrix <- function(var1){
  tmp <- DNAStringSet()
  tmp <- replicate(num_of_reads,generate_reads())
  tmp <- do.call(c,tmp)
  return(tmp)
}
start_time <- proc.time()
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores,type = "PSOCK")
clusterEvalQ(cl,library(Biostrings))
clusterExport(cl, varlist=c("one_subseq","two_subseq","generate_reads","num_of_reads","en1","length_of_reads","name_chr4"))
reads_matrix <- parSapply(cl, 1:num_of_lab,generate_reads_matrix)
#reads_matrix <- parLapply(cl,v,generate_reads_matrix)
stopCluster(cl)
end_time <- proc.time()
print("The time that generating the fragment matrix spends is:")
print(end_time-start_time)
#start_time <- proc.time()
#reads_matrix <- replicate(num_of_lab,generate_reads_matrix())
#end_time <- proc.time()
#end_time-start_time
print("The third step has finished.")
print("You can find one files in the current fold, snp_info.csv.")
rm(one_subseq)
print("The one_subseq has been removed.")
rm(two_subseq)
print("The two_subseq has been removed.")
print("####################################################################")

#print("The fourth step is saving the segment matrix as fastq.")
#get_shortreads <- function(var2){
#  tmp1 <- DNAStringSet(substring(var2[1],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp2 <- DNAStringSet(substring(var2[2],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp3 <- DNAStringSet(substring(var2[3],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp4 <- DNAStringSet(substring(var2[4],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp5 <- DNAStringSet(substring(var2[5],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp6 <- DNAStringSet(substring(var2[6],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp7 <- DNAStringSet(substring(var2[7],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp8 <- DNAStringSet(substring(var2[8],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp9 <- DNAStringSet(substring(var2[9],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp10 <- DNAStringSet(substring(var2[10],seq(1,length_of_reads-127,128),seq(128,length_of_reads,128)))
#  tmp  <- c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10)
#  names(tmp) <- replicate(1600,name_chr4)
#  return(tmp)
#}

#short_matrix <- lapply(reads_matrix, get_shortreads)
#short_matrix <- do.call(c,short_matrix)
#setwd("C:/Users/liutongtong/Desktop/project1/1000x10/test1")
#writeXStringSet(short_matrix, filepath ="segmatrix.fastq", append=FALSE,compress=FALSE, compression_level=NA, format="fastq")
#print("The fourth step has finished.")
#print("####################################################################")
#print("This section has finished, Please enjoy it!")
#print("God bless you!!!!!!!!!!!!!!")


print("The fourth step is generating reads according to pool_method.txt. Maybe this step needs a lot of time.")
pool_method <- read.table("pool_method.txt",sep = ';',header = F)
nrow_of_pool <- dim(pool_method)[1]
ncol_of_pool <- dim(pool_method)[2]
pool_method <- matrix(unlist(pool_method),nrow = nrow_of_pool,ncol = ncol_of_pool )
pool_method <- pool_method[,-ncol_of_pool]
ncol_of_pool <- ncol_of_pool-1
#pool_method <- matrix(1:35,nrow = 7)
#nrow_of_pool <- dim(pool_method)[1]
#ncol_of_pool <- dim(pool_method)[2]
t <- 10
#length_of_getreads <- 5
length_of_getreads <- 250
en3 <- length_of_reads+1-length_of_getreads
num_of_eachreads  <- ceiling(length_of_reads/length_of_getreads)
for (i in 1:3) {
  i <- 1
  print(paste("Th",i,"pool is starting."))
  start_time <- proc.time()
  perpool_experiment <- list()
  for (k in 1:ncol_of_pool) {
    num <- pool_method[i,k]
    perpool_experiment <- c(perpool_experiment,reads_matrix[[num]])
  }
  perpool_experiment <- do.call(c,perpool_experiment)
  len_of_perpool <- length(perpool_experiment)
  getreads <- function(var1){
    tmp_seq <- sample(1:len_of_perpool,size = 1,replace = T)
    tmp_start <- sample(1:en3,size = 1,replace = T)
    tmp <- DNAStringSet(perpool_experiment[tmp_seq], start = tmp_start, end = tmp_start-1+length_of_getreads)
    return(tmp)
  }
  #num_sample <- t*num_of_eachreads*ncol_of_pool
  num_sample <- t*num_of_eachreads*ncol_of_pool*num_of_reads
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores, type = "FORK")
  #cl <- makeCluster(n_cores, type = "PSOCK")
  clusterEvalQ(cl,library(Biostrings))
  clusterExport(cl, varlist=c("perpool_experiment","len_of_perpool","en3","length_of_getreads"))
  perpool_experiment_getreads <- parSapply(cl, 1:num_sample,getreads)
  stopCluster(cl)
  perpool_experiment_getreads <- do.call(c,perpool_experiment_getreads)
  writeXStringSet(perpool_experiment_getreads, filepath =paste("th",i,"pool",sep = "."), append=FALSE,compress=FALSE, compression_level=NA, format="fastq")
  print(paste("Th",i,"pool has finished"))
  end_time <- proc.time()
  print(paste("The time that th",i,"pool spends is:"))
  print(end_time-start_time)
  rm(perpool_experiment_getreads)
  rm(perpool_experiment)
}
print("The fourth step has finished.")
print("####################################################################")
print("This section has finished, Please enjoy it!")
print("God bless you!!!!!!!!!!!!!!")


