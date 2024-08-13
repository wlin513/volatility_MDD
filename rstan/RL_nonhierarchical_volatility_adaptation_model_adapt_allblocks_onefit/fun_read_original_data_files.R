fun_read_original_data_files <- function(file.name, ...){
  lines <- scan(file.name, what="character", sep="\n")
  first.line <- min(grep("trial_number", lines))
  sub_num<-substr(lines[1],13,16)
  blktype<-substr(lines[1],42,43)
  nsubnum=as.numeric(sub_num)
  nblktype=as.numeric(blktype)
  ret_list <- list("data"=read.delim(textConnection(lines), skip=first.line-1),"subnum"=nsubnum,'blktype'=nblktype)
  return(ret_list)
  }
