library(stringr)

setwd("C:/Users/liutongtong/Desktop/project1/1000x10/out1")

data1 <- read.table(file = 'bac_4x.out', sep = ' ', fill = TRUE)
data_f <- as.character(data1[[1]])
data_s <- strsplit(data_f,split = "\t")
num <- length(data_f)
pattern <- "[0-9]+"

for (i in 1:num) {
  #i <- 1
  aloci <- str_extract_all(data_f[i], regex(pattern))
  aloci <- as.numeric(aloci[[1]])
  dis_loci <- aloci[3:length(aloci)]-aloci[2:(length(aloci)-1)]
  dis_ind <- which(dis_loci < 6586)+1
  a <- 1
  b <- 2
  tmp <- c(dis_ind[1])
  ind <- list()
  n <- length(dis_ind)
  while (a<=n) {
    if(a == n && n == 1){
      ind <- c(ind,list(tmp))
      #print("1")
      break
    }else if(a == n && n != 1){
      if(dis_ind[a]!=(dis_ind[a-1]+1)){
        ind <- c(ind, list(dis_ind[a]))
        #print("2")
        break
      }else{
        ind <- c(ind, list(tmp))
        break
      }
    }
    if(dis_ind[b]==(dis_ind[a]+1)){
      tmp <- c(tmp,dis_ind[b])
      a <- a+1
      b <- b+1
      #print("3")
    }else{
      ind <- c(ind,list(tmp))
      a <- b
      b <- a + 1
      #print("4")
      tmp <- c(dis_ind[a])
    }
  }
  
  get_ind <- function(var){
    var1 <- max(var)+1
    tmp <- c(var,var1)
    return(tmp)
  }
  
  f <- lapply(ind, get_ind)
  
  #f_test <- list(f[[16]],f[[39]],f[[34]],237:249)
  
  deal_less_20K <- function(vec) {
    if(length(vec)<=2){
      return(vec)
    }else{
      n2 <- length(vec)
      for (k in n2:3) {
        distance <- aloci[vec[k]]-aloci[vec[1]]
        if(distance > 20000)
          vec <- vec[-k]
      }
      return(vec)
    }
  }
  
  f <- lapply(f, deal_less_20K)
  
  get_loci <- function(var){
    return(data_s[[i]][var])
  }
  
  f1 <- lapply(f,get_loci)
  ind_order <- paste("The ",i," pool")
  f1 <- append(f1, list(ind_order), 0)
  f2 <- lapply(f1, write, "get.fragments1.bac_4x.txt", sep=" ",append=TRUE, ncolumns=1000)
}


