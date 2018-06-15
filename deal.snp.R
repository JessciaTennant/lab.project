library(data.table)
#setwd("C:/Users/liutongtong/Desktop/project1/1000x10/pipeline/the_result")
setwd("C:/Users/liutongtong/Desktop/project1/1000x10/vcf.collection")
filelist <- list.files(pattern = ".*.txt")
datalist <- lapply(filelist, function(x)read.table(x,header = T,col.names = c("loci","orig_base","mutation_base")))
data <- do.call("rbind",datalist)
get_snp <- as.data.table(data)[, list(list(.I)), by = data]
snp_info <- read.table("snp_info.csv",sep = ",",header = T,col.names = c("number","loci","orig_base","mutation_base"))
detect_snp_loci <- intersect(snp_info[[2]],get_snp[[1]])

deal_snp_info <- file("C:/Users/liutongtong/Desktop/project1/1000x10/pipeline/deal_snp_info.txt")
writeLines(c(paste("We has generated ",length(snp_info[[1]])," SNPs totally."),
             paste("We totally detect ",length(data[[1]])," SNPs in 70 pools, including the duplicate values."),
             paste("We get ",length(get_snp[[1]])," different SNPs."),
             paste("We get ",length(detect_snp_loci)," SNPs from we generated."),
             paste("So we get ",length(detect_snp_loci)/length(snp_info[[1]])," SNPs from we generated.")),deal_snp_info)
close(deal_snp_info)

count_perpool_info <- function(var1){
  var1 <- as.data.table(var1)
  var2 <- var1[, .N, by=loci]
  get_number <- match(var2$loci,snp_info$loci)
  var2$loci <- get_number
  return(var2)
}

datalist1 <- lapply(datalist, count_perpool_info)

for (i in 1:length(datalist1)) {
  write.csv(datalist1[[i]], file = paste("C:/Users/liutongtong/Desktop/project1/1000x10/pipeline/count/th_",i,"_pool_count"))
}
