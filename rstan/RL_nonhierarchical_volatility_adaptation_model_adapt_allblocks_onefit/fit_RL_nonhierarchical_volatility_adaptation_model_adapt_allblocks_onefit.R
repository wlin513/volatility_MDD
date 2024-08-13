# script to read in data from 2017_clinical_study single subject data and fit parameters

rm(list=ls())


# configure filepaths -----------------------------------------------------

setwd("D:/2017_clinical_study/code/rstan")
#where the data files are
data_dir<-"D:/2017_clinical_study/data"
dest_dir<-"D:/2017_clinical_study/data/rstan/RL_nonhierarchical_volatility_adaptation_model_adapt_allblocks_onefit"
stan_dir<-"D:/2017_clinical_study/code/rstan/RL_nonhierarchical_volatility_adaptation_model_adapt_allblocks_onefit"
fig_dir<-"D:/2017_clinical_study/tmp_fig/rstan/RL_nonhierarchical_volatility_adaptation_model_adapt_allblocks_onefit"

ifelse(!dir.exists(dest_dir), dir.create(dest_dir), NA)
ifelse(!dir.exists(fig_dir), dir.create(fig_dir),NA)
#some user functions
source(file.path(stan_dir,"fun_read_original_data_files.R"))
source(file.path(stan_dir,"fun_parse_data_into_blocks.R"))
source(file.path(stan_dir,"inv_logit.R"))
source(file.path(stan_dir,"logit.R"))
source(file.path(stan_dir,"read_ests_vol_adpt_model_adapt_allblocks.R"))

#libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# configure task details --------------------------------------------------

#num trials in block
ntrials=80
nabandon=5
nblk=3

#the block order information
blockorders<-matrix(c(1,3,2, 1,2,3),nrow=3,ncol=2)# 1 for both volatile 2 for win volatile 3 for loss volatile
#read in the list of subjects
subjects<-read.table(file=file.path(data_dir,"rsran_sub_session1_data_list_exc_1004_1039_2034.txt"))

#total number of subjects
nsubs=nrow(subjects)


#create arrays for data. For the prescan block, there are 120 trials. For the scan data, there are 80 trials per block and 4 blocks

win_outcome<-array(dim=c(ntrials,nblk,nsubs))
choice<-array(dim=c(ntrials,nblk,nsubs))
loss_outcome<-array(dim=c(ntrials,nblk,nsubs))
includeTrial<-array(dim=c(ntrials,nblk,nsubs))
# loop through subjects getting their data ----------------------------------------
              
for(sn in 1:nsubs)
    {
  sub_data<-fun_parse_data_into_blocks(fun_read_original_data_files(file.path(data_dir,subjects[sn,1]),sep=""))

  #sort data for different block order
  blockorder=blockorders[,sub_data$blktype]
   i=1
   
   for(block in blockorder)
        { print(block)
          print(i)
          win_outcome[,block,sn]<-data.matrix(sub_data$data[paste("block",i,"_wins",sep="")])
          choice[,block,sn]<-data.matrix(sub_data$data[paste("block",i,"_choice",sep="")])
          loss_outcome[,block,sn]<-data.matrix(sub_data$data[paste("block",i,"_loss",sep="")])
          ifmakechoice=is.finite(data.matrix(sub_data$data[paste("block",i,"_RT",sep="")]))
          includet=as.numeric(ifmakechoice)
          includet[1:nabandon]=0
          includeTrial[,block,sn]<-includet
          i=i+1
      }
}
# fit stan ----------------------------------------------------------------



#scan blocks
stan_data<-list(ntr=ntrials,nblk=nblk, nsub=nsubs, opt1Chosen=choice, opt1winout=win_outcome, opt1lossout=loss_outcome,includeTrial=includeTrial)
fit_scandata <- stan(file = file.path(stan_dir,'RL_nonhierarchical_volatility_adaptation_model_adapt_allblocks_onefit.stan'), data = stan_data, iter = 10000, chains = 4, control=list(adapt_delta=0.99, max_treedepth = 15))
save(fit_scandata,file=file.path(dest_dir,"RL_nonhierarchical_volatility_adaptation_model_adapt_allblocks_onefit.RData"))
library("shinystan")
my_sso <- launch_shinystan(fit_scandata)



# figures -----------------------------------------------------------------
#read data and save as mat file
#load(file=file.path(dest_dir,"RL_nonhierarchical_volatility_adaptation_model_adapt_allblocks_onefit.RData"))

IDs=sub('\\/.*',"",subjects$V1)
blocknames=c("both volatile","win volatile","loss volatile")
alldata=read_ests(fit_scandata,IDs,blocknames)

library(ggplot2)
library(plyr)
library(reshape2)
betas=subset(alldata,alldata$variables=="beta")
win_vol_adpt_tr=subset(alldata,alldata$variables=="win_vol_adpt_tr")
loss_vol_adpt_tr=subset(alldata,alldata$variables=="loss_vol_adpt_tr")
winalphas=subset(alldata,alldata$variables=="winalpha")
lossalphas=subset(alldata,alldata$variables=="lossalpha")
allalphas=rbind(winalphas,lossalphas)

lgi_alphas=allalphas
lgi_alphas$values=logit(allalphas$values)


alphameans <- ddply(lgi_alphas, c("block", "variables"), summarise,
               mean=inv_logit(mean(values)),se=inv_logit(mean(values)+sd(values)/sqrt(length(values)))-inv_logit(mean(values)))
palpha<-ggplot(alphameans,aes(x=block, y=mean, fill=variables)) +
            xlab("blocks")+ylab("alphas")+
            geom_bar(position=position_dodge(), stat="identity",
            size=.3) +      # Thinner lines
            scale_fill_manual(values=c("lightgreen","red"))+
            theme_minimal()+
            geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) 

palpha

ggsave("alphas_RL_nonhierarchical_vol_adpt_model_adapt_onefit.eps", plot = palpha, device = "eps", path = fig_dir,
       scale = 1, width = 20, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE)


                          
# save data in a mat file -------------------------------------------------
library(R.matlab)
writeMat(file.path(dest_dir,'Result_vol_adpt_model_adapt.mat'),winalphas=winalphas,lossalphas=lossalphas,betas=betas,win_vol_adpt_tr=win_vol_adpt_tr,loss_vol_adpt_tr=loss_vol_adpt_tr)