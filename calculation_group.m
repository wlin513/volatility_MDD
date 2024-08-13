%% 
clear all;
clc;
%%
cd('/home/wlin/Documents/2017_clinical_study/code');
getfolders;
%%
addpath('plot_code/')
addpath('useful_fun/')
addpath("model-rescorla_wagner_1lr1b/")
addpath("model-rescorla_wagner_2lr1b/")
addpath("model-rescorla_wagner_2lr1b_choice_stickness/")
addpath("model-rescorla_wagner_2lr2b/")
addpath("model-rescorla_wagner_PNPE_2lr1b/")
addpath("model-rescorla_wagner_2lr1b_bias/")
addpath("model-rescorla_wagner_PNPE_2lr1b_opt1/")    
addpath("model-rescorla_wagner_PNPE_2lr2b/")
addpath("model-rescorla_wagner_PNPE_2lr2b_chosen/")
addpath("model-rescorla_wagner_4lr1b_chosen/")
%% load subject data pool

sublist.G1=1001:1041;sublist.G1(sublist.G1==1036)=[];
sublist.G2=2001:2042; sublist.G2(sublist.G2==2033)=[];sublist.G2(sublist.G2==2019)=[];
sublist.G3=3001:3041; sublist.G3(sublist.G3==3033)=[];
% exclusion: exclude participant who chose shape 1 or left side ???
%G1
 sublist.G1(sublist.G1==1004)=[];%who chose shape 1 for 66% 0% 97% for the 3 blocks respectively in session1
 sublist.G1(sublist.G1==1039)=[];%who chose left side for 38.75% 0% 42.5% for the 3 blocks respectively in session1
% sublist.G1(sublist.G1==1017)=[];%who chose shape 1 for 45% 87.5% 1.25% for the 3 blocks respectively in session 1
% sublist.G1(sublist.G1==1010)=[];%who chose shape 1 for 32.5% 85% 95% for the 3 blocks respectively in session 1
% sublist.G1(sublist.G1==1016)=[];%who chose shape 1 for 32.5% 56.25% 93.75% for the 3 blocks respectively in session 1
% sublist.G1(sublist.G1==1021)=[];%who chose shape 1 for 42.5% 93.75% 68.75% for the 3 blocks respectively in session 1

%G2
sublist.G2(sublist.G2==2034)=[];%who chose shape 1 for 0% 93% 85% for the 3 blocks respectively in session1
%sublist.G2(sublist.G2==2032)=[];%who chose shape 1 for 53.75% 3.75% 21.25% for the 3 blocks respectively in session1

%G3
%sublist.G3(sublist.G3==3034)=[];%who chose shape 1 for 23.75% 47.5% 90% for the 3blocks respectively in the session1

sublist.all=[sublist.G1,sublist.G2,sublist.G3];
groups=[ones(size(sublist.G1,2),1);ones(size(sublist.G2,2),1)*2;ones(size(sublist.G3,2),1)*3];
%% get training info
load([datadir,'finaldata.mat'])
for subset={'G1','G2','G3','all'}
   subspool=sublist.(subset{1});
    for ss=1:size(subspool,2)
        subi=find(strcmp(dataout.subject_names,['GB-25-',num2str(subspool(ss))]));
       
        type=extractBefore(dataout.subject_data(subi).IBLT.d1.file,4);
        if strcmp(type,'pos')
             trainingtype.(subset{1})(ss,1)=1;
        else
            if strcmp(type,'con')
                             trainingtype.(subset{1})(ss,1)=2;
            end
        end
    end
end
%%
QIDS=get_p1vital_info(dataout,sublist,trainingtype,'QIDS');
sSTAI=get_p1vital_info(dataout,sublist,trainingtype,'sSTAI');
tSTAI=get_p1vital_info(dataout,sublist,trainingtype,'tSTAI');
SHAPS=get_p1vital_info(dataout,sublist,trainingtype,'SHAPS');
DAQ=get_p1vital_info(dataout,sublist,trainingtype,'DAQ');
RRS=get_p1vital_info(dataout,sublist,trainingtype,'RRS');

%% experiment sdata(i).information setting
tn=80;
blkn=3;
abandontn=10; % abandon first excnum trials when calculating posterior probability
start=[0.5,0.5];
blkname={'both volatile','win volatile','loss volatile'};
blkindex=[1,3,2;1,2,3];%1 for both volatile 2 for win volatile 3 for loss volatile
resp_made=ones(tn,1);
subspool=sublist.all;
%% model fitting for  1lr 1beta model
for session={'session1','session2','session3'}
    %
    for ss=1:size(subspool,2)
    %
    oridata=importdata([datadir,num2str(subspool(ss)),'/',num2str(subspool(ss)),'_',strrep(session{1},'n','n_'),'_vol_train_log.dat']);
    sort_data_into_blocks
        %fit models
        for i=1:blkn
        sdata(i).choice(sdata(i).choice==2)=0;
        tmp_result.(session{1}).all.all(ss,i,:)=fit_linked_1lr_1beta_il(sdata(i).information,sdata(i).choice, start, abandontn);
        end
    end
end
%get full structure 
tmp_result=get_full_structure(tmp_result,trainingtype,sublist);
%get fields for each variable
for session={'session1','session2','session3'}
    for subset={'G1','G2','G3','all'}
        for trainset={'pos','con','all'}
            for i=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),1)
                for j=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),2)
                      alpha_1lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha');
                      beta_1lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_beta'); 
                      BIC_1lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'BIC'); 
                      AIC_1lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'AIC'); 
                end
            end
        end
    end
end
clearvars tmp_result
% %plot bargraph for mean alpha and beta for 1lr 1beta model    
% for session={'session1','session2','session3'}
%     for group={'G1','G2','G3'}
%         plot_bargarph_alphas_rescorla_wagner_2lr_1b
%         %plot_bargarph_betas_rescorla_wagner_2lr_1b
%         alpha_diff_2lr_1b.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_1b.(session{1}).(group{1}).all)-inv_logit(alpha_loss_2lr_1b.(session{1}).(group{1}).all);
%     end
% end
% %%
% plot_bargarph_3group_comparison(alpha_rew_2lr_1b,'alpha rew','session1','learning rate','auto',blkname,figdir)
% plot_bargarph_3group_comparison(alpha_loss_2lr_1b,'alpha loss','session1','learning rate','auto',blkname,figdir)
% plot_bargarph_3group_comparison(alpha_diff_2lr_1b,'alpha diff','session1','learning rate','auto',blkname,figdir)
%% model fitting for  2lr 1beta model
for session={'session1','session2','session3'}
    for ss=1:size(subspool,2)
    oridata=importdata([datadir,num2str(subspool(ss)),'/',num2str(subspool(ss)),'_',strrep(session{1},'n','n_'),'_vol_train_log.dat']);
    sort_data_into_blocks
        %fit models
        for i=1:blkn
        sdata(i).choice(sdata(i).choice==2)=0;
        %2alphas1beta
        tmp_result.(session{1}).all.all(ss,i,:)=fit_linked_2lr_1beta_il(sdata(i).information,sdata(i).choice, start, abandontn);
        end
    end
end
%get full structure 
tmp_result=get_full_structure(tmp_result,trainingtype,sublist);
%get fields for each variable
for session={'session1','session2','session3'}
    for subset={'G1','G2','G3','all'}
        for trainset={'pos','con','all'}
            for i=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),1)
                for j=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),2)
                      alpha_rew_2lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha_rew');
                      alpha_loss_2lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha_loss');
                      beta_2lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_beta'); 
                      BIC_2lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'BIC'); 
                      AIC_2lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'AIC'); 
                end
            end
        end
    end
end
clearvars tmp_result
%%plot bargraph for mean alpha and beta for 1lr 1beta model    
for session={'session1','session2','session3'}
    for group={'G1','G2','G3'}
        plot_2lrs(alpha_rew_2lr_1b,alpha_loss_2lr_1b,group{1},session{1},'2lr_1b','auto')
        %plot_bargarph_betas_rescorla_wagner_2lr_1b
        alpha_diff_2lr_1b.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_1b.(session{1}).(group{1}).all)-inv_logit(alpha_loss_2lr_1b.(session{1}).(group{1}).all);
        alpha_rew_bothv_lossv_2lr_1b.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_1b.(session{1}).(group{1}).all(:,1))-inv_logit(alpha_rew_2lr_1b.(session{1}).(group{1}).all(:,3));
        alpha_loss_bothv_winv_2lr_1b.(session{1}).(group{1}).all=inv_logit(alpha_loss_2lr_1b.(session{1}).(group{1}).all(:,1))-inv_logit(alpha_loss_2lr_1b.(session{1}).(group{1}).all(:,2));
        alpha_rew_winv_lossv_2lr_1b.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_1b.(session{1}).(group{1}).all(:,2))-inv_logit(alpha_rew_2lr_1b.(session{1}).(group{1}).all(:,3));
        alpha_loss_lossv_winv_2lr_1b.(session{1}).(group{1}).all=inv_logit(alpha_loss_2lr_1b.(session{1}).(group{1}).all(:,3))-inv_logit(alpha_loss_2lr_1b.(session{1}).(group{1}).all(:,2));

    end
end
%
plot_bargarph_3group_comparison_lrs(alpha_rew_2lr_1b,'alpha rew','session1','2lr_1b','auto','northeast',blkname,figdir)
plot_bargarph_3group_comparison_lrs(alpha_loss_2lr_1b,'alpha loss','session1','2lr_1b','auto','northwest',blkname,figdir)
plot_bargarph_3group_comparison(alpha_diff_2lr_1b,'alpha diff','session1','learning rate','auto',blkname,figdir)
plot_bargarph_3group_comparison_betas(beta_2lr_1b,'beta','session1','2lr_1b','auto','northeast',blkname,figdir)
plot_bargarph_3group_comparison_mean(alpha_rew_bothv_lossv_2lr_1b,'session1','win alpha bothv vs lossv','learning rate','auto',figdir)
plot_bargarph_3group_comparison_mean(alpha_loss_bothv_winv_2lr_1b,'session1','loss alpha bothv vs winv','learning rate','auto',figdir)
plot_bargarph_3group_comparison_mean(alpha_rew_winv_lossv_2lr_1b,'session1','win alpha winv vs lossv','learning rate','auto',figdir)
plot_bargarph_3group_comparison_mean(alpha_loss_lossv_winv_2lr_1b,'session1','loss alpha lossv vs winv','learning rate','auto',figdir)

%% model fitting for  2lr 1beta option bias model
for session={'session1','session2','session3'}
    for ss=1:size(subspool,2)
    oridata=importdata([datadir,num2str(subspool(ss)),'/',num2str(subspool(ss)),'_',strrep(session{1},'n','n_'),'_vol_train_log.dat']);
    sort_data_into_blocks
        %fit models
        for i=1:blkn
        sdata(i).choice(sdata(i).choice==2)=0;
        %2alphas1beta+1bias term
        tmp_result.(session{1}).all.all(ss,i,:)=fit_linked_2lr_beta_add(sdata(i).information,sdata(i).choice, start, abandontn);
        end
    end
end
%get full structure 
tmp_result=get_full_structure(tmp_result,trainingtype,sublist);
%get fields for each variable
for session={'session1','session2','session3'}
    for subset={'G1','G2','G3','all'}
        for trainset={'pos','con','all'}
            for i=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),1)
                for j=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),2)
                      alpha_rew_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha_rew');
                      alpha_loss_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha_loss');
                      beta_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_beta'); 
                      BIC_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'BIC'); 
                      AIC_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'AIC'); 
                end
            end
        end
    end
end
clearvars tmp_result
%%plot bargraph for mean alpha and beta for 1lr 1beta model    
for session={'session1','session2','session3'}
    for group={'G1','G2','G3','all'}
        plot_2lrs(alpha_rew_2lr_1b_bias,alpha_loss_2lr_1b_bias,group{1},session{1},'2lr_1b_bias','auto')
        alpha_diff_2lr_1b_bias.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_1b_bias.(session{1}).(group{1}).all)-inv_logit(alpha_loss_2lr_1b_bias.(session{1}).(group{1}).all);
        alpha_rew_bothv_lossv_2lr_1b_bias.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_1b_bias.(session{1}).(group{1}).all(:,1))-inv_logit(alpha_rew_2lr_1b_bias.(session{1}).(group{1}).all(:,3));
        alpha_loss_bothv_winv_2lr_1b_bias.(session{1}).(group{1}).all=inv_logit(alpha_loss_2lr_1b_bias.(session{1}).(group{1}).all(:,1))-inv_logit(alpha_loss_2lr_1b_bias.(session{1}).(group{1}).all(:,2));
        alpha_rew_winv_lossv_2lr_1b_bias.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_1b_bias.(session{1}).(group{1}).all(:,2))-inv_logit(alpha_rew_2lr_1b_bias.(session{1}).(group{1}).all(:,3));
        alpha_loss_lossv_winv_2lr_1b_bias.(session{1}).(group{1}).all=inv_logit(alpha_loss_2lr_1b_bias.(session{1}).(group{1}).all(:,3))-inv_logit(alpha_loss_2lr_1b_bias.(session{1}).(group{1}).all(:,2));

    end
end
%
plot_bargarph_3group_comparison_lrs(alpha_rew_2lr_1b_bias,'alpha rew','session1','2lr_1b_bias','auto','northeast',blkname,figdir)
plot_bargarph_3group_comparison_lrs(alpha_loss_2lr_1b_bias,'alpha loss','session1','2lr_1b_bias','auto','northwest',blkname,figdir)
plot_bargarph_3group_comparison(alpha_diff_2lr_1b_bias,'alpha diff','session1','learning rate','auto',blkname,figdir)
plot_bargarph_3group_comparison_betas(beta_2lr_1b_bias,'beta','session1','2lr_1b_bias','auto','northeast',blkname,figdir)
plot_bargarph_3group_comparison_mean(alpha_rew_bothv_lossv_2lr_1b_bias,'session1','win alpha bothv vs lossv','learning rate','auto',figdir)
plot_bargarph_3group_comparison_mean(alpha_loss_bothv_winv_2lr_1b_bias,'session1','loss alpha bothv vs winv','learning rate','auto',figdir)
plot_bargarph_3group_comparison_mean(alpha_rew_winv_lossv_2lr_1b_bias,'session1','win alpha winv vs lossv','learning rate','auto',figdir)
plot_bargarph_3group_comparison_mean(alpha_loss_lossv_winv_2lr_1b_bias,'session1','loss alpha lossv vs winv','learning rate','auto',figdir)

anova1(alpha_loss_lossv_winv_2lr_1b_bias.session1.all.all,groups)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),sSTAI.session1.all.all,'mean beta','state anxiety',...
    'scatterplot_mean_beta_corr_sSTAI','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),tSTAI.session1.all.all,'mean beta','trait anxiety',...
    'scatterplot_mean_beta_corr_tSTAI','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),QIDS.session1.all.all,'mean beta','QIDS',...
    'scatterplot_mean_beta_corr_QIDS','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),SHAPS.session1.all.all,'mean beta','SHAPS',...
    'scatterplot_mean_beta_corr_SHAPS','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),DAQ.session1.all.all,'mean beta','DAQ',...
    'scatterplot_mean_beta_corr_DAQ','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),RRS.session1.all.all,'mean beta','RRS',...
    'scatterplot_mean_beta_corr_RRS','Pearson',groups,figdir)
blk=1;
[r,p]=partialcorr(inv_logit(alpha_rew_2lr_1b_bias.session1.all.all(:,blk)),sSTAI.session1.all.all,mean(log(beta_2lr_1b_bias.session1.all.all),2))
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),mean(switchprob_win.session1.all.all,2),'mean beta','win switchprob',...
    'scatterplot_mean_beta_corr_win_switchprob','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),mean(switchprob_nowin.session1.all.all,2),'mean beta','nowin switchprob',...
    'scatterplot_mean_beta_corr_nowin_switchprob','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),mean(switchprob_noloss.session1.all.all,2),'mean beta','noloss switchprob',...
    'scatterplot_mean_beta_corr_noloss_switchprob','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),mean(switchprob_loss.session1.all.all,2),'mean beta','loss switchprob',...
    'scatterplot_mean_beta_corr_loss_switchprob','Pearson',groups,figdir)

scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),mean(inv_logit(alpha_rew_2lr_1b_bias.session1.all.all),2),...
    'mean beta','mean win alpha','scatterplot_mean_beta_corr_alpha_rew_2lr_1b_bias','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_2lr_1b_bias.session1.all.all),2),mean(inv_logit(alpha_loss_2lr_1b_bias.session1.all.all),2),...
    'mean beta','mean loss alpha','scatterplot_mean_beta_corr_alpha_loss_2lr_1b_bias','Pearson',groups,figdir)
scatterplot_3groups(mean(inv_logit(alpha_rew_2lr_1b_bias.session1.all.all),2),mean(inv_logit(alpha_loss_2lr_1b_bias.session1.all.all),2),...
    'mean win alpha','mean loss alpha','scatterplot_win_alpha_corr_alpha_loss_2lr_1b_bias','Pearson',groups,figdir)

blk=3;
scatterplot_3groups(log(beta_2lr_1b_bias.session1.all.all(:,blk)),inv_logit(alpha_rew_2lr_1b_bias.session1.all.all(:,blk)),...
    'beta',' win alpha','scatterplot_mean_beta_corr_alpha_rew_2lr_1b_bias','Pearson',groups,figdir)
scatterplot_3groups(log(beta_2lr_1b_bias.session1.all.all(:,blk)),inv_logit(alpha_loss_2lr_1b_bias.session1.all.all(:,blk)),...
    'beta','loss alpha','scatterplot_mean_beta_corr_alpha_loss_2lr_1b_bias','Pearson',groups,figdir)
scatterplot_3groups(inv_logit(alpha_rew_2lr_1b_bias.session1.all.all(:,blk)),inv_logit(alpha_loss_2lr_1b_bias.session1.all.all(:,blk)),...
    'win alpha','loss alpha','scatterplot_win_alpha_corr_alpha_loss_2lr_1b_bias','Pearson',groups,figdir)

scatterplot_3groups(mean(inv_logit(alpha_rew_2lr_1b.session1.all.all),2),mean(inv_logit(alpha_loss_2lr_1b.session1.all.all),2),...
    'mean win alpha','mean loss alpha','-','Pearson',groups,figdir)
scatterplot_3groups(mean(log(beta_loss_2lr_2b.session1.all.all),2),mean(inv_logit(alpha_rew_2lr_2b.session1.all.all),2),...
    'mean beta','mean loss alpha','-','Pearson',groups,figdir)
% scatterplot_3groups(log(beta_2lr_1b_bias.session1.all.all(:,1)),sSTAI.session1.all.all,'mean beta','state anxiety',...
%     'scatterplot_bothv_beta_corr_sSTAI','Pearson',groups,figdir)
%% model fitting for  2lr 2beta model
for session={'session1','session2','session3'}
    for ss=1:size(subspool,2)
    oridata=importdata([datadir,num2str(subspool(ss)),'/',num2str(subspool(ss)),'_',strrep(session{1},'n','n_'),'_vol_train_log.dat']);
    sort_data_into_blocks
        %fit models
        for i=1:blkn
        sdata(i).choice(sdata(i).choice==2)=0;
        %2alphas1beta
        tmp_result.(session{1}).all.all(ss,i,:)=fit_linked_2lr_2beta_il_rev(sdata(i).information,sdata(i).choice, start, abandontn);
        end
    end
end
%get full structure 
tmp_result=get_full_structure(tmp_result,trainingtype,sublist);
%get fields for each variable
for session={'session1','session2','session3'}
    for subset={'G1','G2','G3','all'}
        for trainset={'pos','con','all'}
            for i=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),1)
                for j=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),2)
                      alpha_rew_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha_rew');
                      alpha_loss_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha_loss');
                      beta_rew_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_beta_rew'); 
                      beta_loss_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_beta_loss'); 
                      BIC_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'BIC'); 
                      AIC_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'AIC'); 
                end
            end
        end
    end
end
clearvars tmp_result
%%plot bargraph for mean alpha and beta for 1lr 1beta model    
for session={'session1','session2','session3'}
    for group={'G1','G2','G3'}
        plot_2lrs(alpha_rew_2lr_2b,alpha_loss_2lr_2b,group{1},session{1},'2lr_2b','auto')
        plot_2betas(beta_rew_2lr_2b,beta_loss_2lr_2b,group{1},session{1},'2lr_2b','auto')
        alpha_diff_2lr_2b.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_2b.(session{1}).(group{1}).all)-inv_logit(alpha_loss_2lr_2b.(session{1}).(group{1}).all);
    end
end
%
plot_bargarph_3group_comparison_lrs(alpha_rew_2lr_2b,'alpha rew','session1','2lr_2b','auto','northeast',blkname,figdir)
plot_bargarph_3group_comparison_lrs(alpha_loss_2lr_2b,'alpha loss','session1','2lr_2b','auto','northwest',blkname,figdir)
plot_bargarph_3group_comparison(alpha_diff_2lr_2b,'alpha diff','session1','learning rate','auto',blkname,figdir)
%% model fitting for  pos neg PE 2lr 2beta model
for session={'session1','session2','session3'}
    for ss=1:size(subspool,2)
    oridata=importdata([datadir,num2str(subspool(ss)),'/',num2str(subspool(ss)),'_',strrep(session{1},'n','n_'),'_vol_train_log.dat']);
    sort_data_into_blocks
        %fit models
        for i=1:blkn
        sdata(i).choice(sdata(i).choice==2)=0;
        %%model with seperate alphas for positive and negative prediction errors
        tmp_result.(session{1}).all.all(ss,i,:)=fit_linked_PNPE_2lr_2beta_chosen(sdata(i).information,sdata(i).choice, start, abandontn);
        end
    end
end
%get full structure 
tmp_result=get_full_structure(tmp_result,trainingtype,sublist);
%get fields for each variable
for session={'session1','session2','session3'}
    for subset={'G1','G2','G3','all'}
        for trainset={'pos','con','all'}
            for i=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),1)
                for j=1:size(tmp_result.(session{1}).(subset{1}).(trainset{1}),2)
                      alpha_PPE_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha_pos');
                      alpha_NPE_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_alpha_neg');
                      beta_rew_PNPE_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_beta_rew'); 
                      beta_loss_PNPE_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'mean_beta_loss'); 
                      BIC_PNPE_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'BIC'); 
                      AIC_PNPE_2lr_2b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(tmp_result,session{1},subset{1},trainset{1},{i,j},'AIC'); 
                end
            end
        end
    end
end
clearvars tmp_result
%%plot bargraph for mean alpha and beta for 1lr 1beta model    
for session={'session1','session2','session3'}
    for group={'G1','G2','G3'}
        plot_2lrs(alpha_PPE_2lr_2b,alpha_NPE_2lr_2b,group{1},session{1},'PNPE_2lr_2b','auto')
       plot_2betas(beta_rew_PNPE_2lr_2b,beta_loss_PNPE_2lr_2b,group{1},session{1},'PNPE_2lr_2b','auto')
        %alpha_diff_2lr_2b.(session{1}).(group{1}).all=inv_logit(alpha_rew_2lr_2b.(session{1}).(group{1}).all)-inv_logit(alpha_loss_2lr_2b.(session{1}).(group{1}).all);
    end
end
%
plot_bargarph_3group_comparison_lrs(alpha_PPE_2lr_2b,'alpha pos PE','session1','PNPE_2lr_2b','auto','northeast',blkname,figdir)
plot_bargarph_3group_comparison_lrs(alpha_NPE_2lr_2b,'alpha neg PE','session1','PNPE_2lr_2b','auto','northwest',blkname,figdir)
%plot_bargarph_3group_comparison(alpha_diff_2lr_2b,'alpha diff','session1','learning rate','auto',blkname,figdir)
plot_bargarph_3group_comparison_betas(beta_rew_PNPE_2lr_2b,'rew beta','session1','PNPE_2lr_1b','auto','northeast',blkname,figdir)
plot_bargarph_3group_comparison_betas(beta_loss_PNPE_2lr_2b,'loss beta','session1','PNPE_2lr_1b','auto','northeast',blkname,figdir)

    %    result_PNPE_2lr_2b.(session{1}).all.all(ss,i,:)=fit_linked_PNPE_2lr_2beta_chosen(sdata(i).information,sdata(i).choice, start, abandontn);
  
%%
save('model_results.mat')
%load('model_results.mat')

%% %stats
T_s=array2table(sublist.all','VariableNames',{'subject'});
T_group=array2table(groups,'VariableNames',{'group'});
T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});
T_alphas_inv=array2table([inv_logit(alpha_rew_2lr_2b.session1.all.all),inv_logit(alpha_loss_2lr_2b.session1.all.all)],...
    'VariableNames',{'rew_alpha_bothv','rew_alpha_winv','rew_alpha_lossv',...
    'loss_alpha_bothv','loss_alpha_winv','loss_alpha_lossv'});
T=[T_s,T_group,T_treatment,T_alphas_inv];
% within=table(repelem({'win';'loss'},3,1),repmat({'winv';'winv';'wins'},2,1),repmat({'lossv';'losss';'lossv'},2,1),...
%      'VariableNames',{'valence','winvol','lossvol'});
within=table(repelem({'win';'loss'},3,1),repmat({'1';'2';'3'},2,1),...
     'VariableNames',{'valence','blktype'});
rm=fitrm(T,'rew_alpha_bothv-loss_alpha_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','valence*blktype');

%
within=table({'1';'2';'3'},...
     'VariableNames',{'blktype'});
rm=fitrm(T,'rew_alpha_bothv-rew_alpha_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','blktype');
%
within=table({'1';'2';'3'},...
     'VariableNames',{'blktype'});
rm=fitrm(T,'loss_alpha_bothv-loss_alpha_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','blktype');
 %oneway anova
 %win
 blk=2;
 [p,tbl,stats]=kruskalwallis([inv_logit(alpha_rew_2lr_2b.session1.G1.all(:,blk));...
     inv_logit(alpha_rew_2lr_2b.session1.G2.all(:,blk));...
     inv_logit(alpha_rew_2lr_2b.session1.G3.all(:,blk))],...
     groups);
 
 %post-hoc ttest
 blk=2;
 [h,p,ci,stats]=ttest2(inv_logit(alpha_rew_2lr_2b.session1.G3.all(:,blk)),inv_logit(alpha_rew_2lr_2b.session1.G2.all(:,blk)))
 
 %loss
 blk=3;
 [p,tbl,stats]=kruskalwallis([inv_logit(alpha_loss_2lr_2b.session1.G1.all(:,blk));...
     inv_logit(alpha_loss_2lr_2b.session1.G2.all(:,blk));...
     inv_logit(alpha_loss_2lr_2b.session1.G3.all(:,blk))],...
     groups);
 
  blk=3;
 [p,tbl,stats]=anova1([pct_loss_alpha.session1.G1.all(:,blk);...
     pct_loss_alpha.session1.G2.all(:,blk);...
     pct_loss_alpha.session1.G3.all(:,blk)],...
     group);
  %post-hoc ttest
 blk=3;
 [h,p,ci,stats]=ttest2(pct_loss_alpha.session1.G1.all(:,blk),pct_loss_alpha.session1.G3.all(:,blk))
%
within=table({'1';'2';'3'},...
     'VariableNames',{'blktype'});
rm=fitrm(T,'winvslosschosen_bothv-winvslosschosen_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','blktype');
%%
for group={'G1','G2','G3','all'}
delta_win.session1.(group{1}).all=inv_logit(alpha_rew_2lr_1b.session1.(group{1}).all(:,2))-inv_logit(alpha_rew_2lr_1b.session1.(group{1}).all(:,3));
delta_loss.session1.(group{1}).all=inv_logit(alpha_loss_2lr_1b.session1.(group{1}).all(:,3))-inv_logit(alpha_loss_2lr_1b.session1.(group{1}).all(:,2));
end
plot_bargarph_3group_comparison_mean(delta_win,'session1','delta-win','delta learning rate','auto',figdir)
plot_bargarph_3group_comparison_mean(delta_loss,'session1','delta-loss','delta learning rate','auto',figdir)
kruskalwallis([delta_win.session1.G1.all;delta_win.session1.G2.all;delta_win.session1.G3.all],groups)
anova1([delta_win.session1.G1.all;delta_win.session1.G2.all;delta_win.session1.G3.all],groups)
[p,h]=ranksum(delta_win.session1.G1.all,delta_win.session1.G3.all)
scatterplot_3groups(SHAPS.session1.all.all,delta_win.session1.all.all,...
    'trait anxiety','relative win learning rates','scatterplot_tSTAI_delta_win_alpha_2lr1b','Pearson',groups,figdir)
[r,p]=partialcorr(tSTAI.session1.all.all,delta_win.session1.all.all,mean(log(beta_2lr_1b.session1.all.all),2))
scatterplot_3groups(tSTAI.session1.all.all,mean(log(beta_2lr_1b.session1.all.all),2),...
    'trait anxiety','mean log beta','scatterplot_tSTAI_mean_beta_2lr1b','Pearson',groups,figdir)
%[r,p]=partialcorr(tSTAI.session1.all.all,mean(log(beta_2lr_1b.session1.all.all),2),delta_win.session1.all.all,)

blk=3;
anova1([mean(alpha_diff_2lr_2b.session1.G1.all,2);mean(alpha_diff_2lr_2b.session1.G2.all,2);mean(alpha_diff_2lr_2b.session1.G3.all,2)],groups)
[h,p]=ttest2(mean(alpha_diff_2lr_2b.session1.G1.all,2),mean(alpha_diff_2lr_2b.session1.G2.all,2))
blk=2;
[h,p]=ttest2(alpha_diff_2lr_2b.session1.G1.all(:,blk),alpha_diff_2lr_2b.session1.G2.all(:,blk))
[h,p]=ttest2(delta_win.session1.G1.all,delta_win.session1.G3.all)
kruskalwallis([delta_win.session1.G1.all,delta_win.session1.G2.all,delta_win.session1.G3.all])
kruskalwallis([delta_loss.session1.G1.all,delta_loss.session1.G2.all,delta_loss.session1.G3.all])
%%

%% sorting and plotting results from 2lr1b bias model
%get full structure 
result_2lr_1b_bias=get_full_structure(result_2lr_1b_bias,trainingtype,sublist);
%get fields for each variable
for session={'session1','session2','session3'}
    for subset={'G1','G2','G3','all'}
        for trainset={'pos','con','all'}
            for i=1:size(result_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1}),1)
                for j=1:size(result_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1}),2)
                      alpha_rew_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(result_2lr_1b_bias,session{1},subset{1},trainset{1},{i,j},'mean_alpha_rew');
                      alpha_loss_2lr_1b_bias.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(result_2lr_1b_bias,session{1},subset{1},trainset{1},{i,j},'mean_alpha_loss');
                      beta_2lr_1b.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(result_2lr_1b_bias,session{1},subset{1},trainset{1},{i,j},'mean_beta');                   
                end
            end
        end
    end
end
clearvars result_2lr_1b_bias
% plot bargraph for mean alpha and beta for 2lr 1beta model    
for session={'session1','session2','session3'}
    for group={'G1','G2','G3'}
        %plot_bargarph_alphas_rescorla_wagner_2lr_1b
        %plot_bargarph_betas_rescorla_wagner_2lr_1b
        alpha_diff_2lr_1b_bias.(session{1}).(group{1}).all=alpha_rew_2lr_1b_bias.(session{1}).(group{1}).all-alpha_loss_2lr_1b.(session{1}).(group{1}).all;
    end
end
%% plot BICs and AICs for each block
modelnames={'2lr2b','2lr1b','1lr1b','2lr1b+bias','PNPE-2lr2b'};
plot_bargarph_AICorBICs('BIC','G3',modelnames,blkname,[50,90],figdir,BIC_2lr_2b,BIC_2lr_1b,BIC_1lr_1b,BIC_2lr_1b_bias,BIC_PNPE_2lr_2b)
plot_bargarph_AICorBICs('AIC','all',modelnames,blkname,[50,80],figdir,AIC_2lr_2b,AIC_2lr_1b,AIC_1lr_1b,AIC_2lr_1b_bias,AIC_PNPE_2lr_2b)
modelnames={'2lr2b','2lr1b','1lr1b'};
plot_bargarph_AICorBICs('BIC','G3',modelnames,blkname,[50,90],figdir,BIC_2lr_2b,BIC_2lr_1b,BIC_1lr_1b)

%% correlate parameters with behavior performance
%load('model_free_results.mat')
blk=3;
scatterplot_3groups(log(beta_2lr_1b_bias.session1.all.all(:,blk)),pct_winvsloss_chosen.session1.all.all(:,blk),...
    'beta','pct-winvsloss-chosen','scatterplot_beta_2lr_1b_bias_corr_pct_winvsloss_chosen','Pearson',groups,figdir)
scatterplot_3groups(log(beta_2lr_1b_bias.session1.all.all(:,blk)),pct_winchosen.session1.all.all(:,blk),...
    'beta','pct-win-chosen','scatterplot_beta_2lr_1b_bias_corr_pct_win_chosen','Pearson',groups,figdir)
scatterplot_3groups(log(beta_2lr_1b_bias.session1.all.all(:,blk)),pct_lossnotchosen.session1.all.all(:,blk),...
    'beta','pct-loss-not-chosen','scatterplot_beta_2lr_1b_bias_corr_pct_loss-not_chosen','Pearson',groups,figdir)

blk=3;
scatterplot_3groups(alpha_rew_2lr_1b_bias.session1.all.all(:,blk),pct_winvsloss_chosen.session1.all.all(:,blk),...
    'win alpha','pct-winvsloss-chosen','scatterplot_win_alpha_2lr_1b_bias_corr_pct_winvsloss_chosen','Spearman',groups,figdir)
scatterplot_3groups(alpha_rew_2lr_1b_bias.session1.all.all(:,blk),pct_winchosen.session1.all.all(:,blk),...
    'win alpha','pct-win-chosen','scatterplot_win_alpha_2lr_1b_bias_corr_pct_win_chosen','Spearman',groups,figdir)
scatterplot_3groups(alpha_rew_2lr_1b_bias.session1.all.all(:,blk),pct_lossnotchosen.session1.all.all(:,blk),...
    'win alpha','pct-loss-not-chosen','scatterplot_win_alpha_2lr_1b_bias_corr_pct_loss-not_chosen','Spearman',groups,figdir)

blk=3;
scatterplot_3groups(alpha_loss_2lr_1b_bias.session1.all.all(:,blk),pct_winvsloss_chosen.session1.all.all(:,blk),...
    'loss alpha','pct-win-chosen','scatterplot_loss_alpha_2lr_1b_bias_corr_pct_winvsloss_chosen','Spearman',groups,figdir)
scatterplot_3groups(alpha_loss_2lr_1b_bias.session1.all.all(:,blk),pct_winchosen.session1.all.all(:,blk),...
    'loss alpha','pct-win-chosen','scatterplot_loss_alpha_2lr_1b_bias_corr_pct_win_chosen','Spearman',groups,figdir)
scatterplot_3groups(alpha_loss_2lr_1b_bias.session1.all.all(:,blk),pct_lossnotchosen.session1.all.all(:,blk),...
    'loss alpha','pct-loss-not-chosen','scatterplot_loss_alpha_2lr_1b_bias_corr_pct_loss_not_chosen','Spearman',groups,figdir)

blk=3;
scatterplot_3groups(alpha_rew_2lr_1b_bias.session1.all.all(:,blk),switchprob_nowin.session1.all.all(:,blk),...
    'win alpha','nowin-switch probability','scatterplot_win_alpha_2lr_1b_bias_corr_pct_nowin_switch','Spearman',groups,figdir)
scatterplot_3groups(alpha_rew_2lr_1b_bias.session1.all.all(:,blk),switchprob_win.session1.all.all(:,blk),...
    'win alpha','win-switch probability','scatterplot_win_alpha_2lr_1b_bias_corr_pct_win_switch','Spearman',groups,figdir)
scatterplot_3groups(alpha_rew_2lr_1b_bias.session1.all.all(:,blk),switchprob_noloss.session1.all.all(:,blk),...
    'win alpha','noloss-switch probability','scatterplot_win_alpha_2lr_1b_bias_corr_pct_noloss_switch','Spearman',groups,figdir)
scatterplot_3groups(alpha_rew_2lr_1b_bias.session1.all.all(:,blk),switchprob_loss.session1.all.all(:,blk),...
    'win alpha','loss-switch probability','scatterplot_win_alpha_2lr_1b_bias_corr_pct_loss_switch','Spearman',groups,figdir)
scatterplot_3groups(alpha_loss_2lr_1b_bias.session1.all.all(:,blk),switchprob_noloss.session1.all.all(:,blk),...
    'loss alpha','noloss-switch probability','scatterplot_loss_alpha_2lr_1b_bias_corr_pct_noloss_switch','Spearman',groups,figdir)
scatterplot_3groups(alpha_loss_2lr_1b_bias.session1.all.all(:,blk),switchprob_loss.session1.all.all(:,blk),...
    'loss alpha','loss-switch probability','scatterplot_loss_alpha_2lr_1b_bias_corr_pct_loss_switch','Spearman',groups,figdir)
scatterplot_3groups(alpha_loss_2lr_1b_bias.session1.all.all(:,blk),switchprob_nowin.session1.all.all(:,blk),...
    'loss alpha','nowin-switch probability','scatterplot_loss_alpha_2lr_1b_bias_corr_pct_nowin_switch','Spearman',groups,figdir)
scatterplot_3groups(alpha_loss_2lr_1b_bias.session1.all.all(:,blk),switchprob_win.session1.all.all(:,blk),...
    'loss alpha','win-switch probability','scatterplot_loss_alpha_2lr_1b_bias_corr_pct_win_switch','Spearman',groups,figdir)

blk=3;
scatterplot_3groups(switchprob_nowin.session1.all.all(:,blk),pct_lossnotchosen.session1.all.all(:,blk),...
    'nowin - switch probability','pct-loss-not-chosen','scatterplot_pct_nowin_switch_corr_pct_loss_not_chosen','Spearman',groups,figdir)
scatterplot_3groups(switchprob_loss.session1.all.all(:,blk),pct_winchosen.session1.all.all(:,blk),...
    'loss - switch probabilibity ','pct-win-chosen','scatterplot_pct_loss_switch_corr_pct_win_chosen','Spearman',groups,figdir)

blk=3;
scatterplot_3groups(alpha_loss_2lr_1b_bias.session1.all.all(:,blk),alpha_rew_2lr_1b_bias.session1.all.all(:,blk),...
    '-','-','scatterplot_pct_win_chosen_2lr_1b_bias_corr_pct_loss_not_chosen','Spearman',groups,figdir)
scatterplot_3groups(switchprob_nowin.session1.all.all(:,blk),pct_lossnotchosen.session1.all.all(:,blk),...
    '-','-','scatterplot_pct_win_chosen_2lr_1b_bias_corr_pct_loss_not_chosen','Spearman',groups,figdir)

blk=3;
scatterplot_3groups(inv_logit(alpha_loss_2lr_1b.session1.all.all(:,blk)),log(beta_2lr_1b.session1.all.all(:,blk)),...
    'loss alpha','beta','scatterplot_loss_alpha_beta_2lr1b_lv','Pearson',groups,figdir)
scatterplot_3groups(inv_logit(alpha_rew_2lr_1b.session1.all.all(:,blk)),log(beta_2lr_1b.session1.all.all(:,blk)),...
    'win alpha','beta','scatterplot_win_alpha_beta_2lr1b_lv','Pearson',groups,figdir)
%%
T_s=array2table(sublist.all','VariableNames',{'subject'});
T_group=array2table(groups,'VariableNames',{'group'});
T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});
T_alphas=array2table([inv_logit(alpha_rew_2lr_1b.session1.all.all(:,2:3)),inv_logit(alpha_loss_2lr_1b.session1.all.all(:,2:3))],...
    'VariableNames',{'win_alpha_winv','win_alpha_lossv',...
                                'loss_alpha_winv','loss_alpha_lossv'});
T=[T_s,T_group,T_treatment,T_alphas];
% within=table(repelem({'win';'loss'},3,1),repmat({'winv';'winv';'wins'},2,1),repmat({'lossv';'losss';'lossv'},2,1),...
%      'VariableNames',{'valence','winvol','lossvol'});
within=table(repelem({'win';'loss'},2,1),repmat({'winv';'lossv'},2,1),...
     'VariableNames',{'valence','blktype'});
rm=fitrm(T,'win_alpha_winv-loss_alpha_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','valence*blktype');
%%
T_s=array2table(sublist.all','VariableNames',{'subject'});
T_group=array2table(groups,'VariableNames',{'group'});
T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});
T_alphas=array2table([inv_logit(alpha_rew_2lr_1b.session1.all.all(:,2:3))],...
    'VariableNames',{'win_alpha_winv','win_alpha_lossv'});
T=[T_s,T_group,T_treatment,T_alphas];
% within=table(repelem({'win';'loss'},3,1),repmat({'winv';'winv';'wins'},2,1),repmat({'lossv';'losss';'lossv'},2,1),...
%      'VariableNames',{'valence','winvol','lossvol'});
within=table(repmat({'winv';'lossv'},1,1),...
     'VariableNames',{'blktype'});
rm=fitrm(T,'win_alpha_winv-win_alpha_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','blktype');
%%
T_s=array2table(sublist.all','VariableNames',{'subject'});
T_group=array2table(groups,'VariableNames',{'group'});
T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});
T_alphas=array2table([inv_logit(alpha_loss_2lr_1b.session1.all.all(:,2:3))],...
    'VariableNames',{'loss_alpha_winv','loss_alpha_lossv'});
T=[T_s,T_group,T_treatment,T_alphas];
% within=table(repelem({'win';'loss'},3,1),repmat({'winv';'winv';'wins'},2,1),repmat({'lossv';'losss';'lossv'},2,1),...
%      'VariableNames',{'valence','winvol','lossvol'});
within=table(repmat({'winv';'lossv'},1,1),...
     'VariableNames',{'blktype'});
rm=fitrm(T,'loss_alpha_winv-loss_alpha_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','blktype');
%% save results to spss
T_s=array2table(sublist.all','VariableNames',{'subject'});
T_group=array2table(groups,'VariableNames',{'group'});
T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});
T_performance=array2table([pct_winchosen.session1.all.all,pct_lossnotchosen.session1.all.all,pct_winvsloss.session1.all.all],...
    'VariableNames',{'winchosen_bothv','winchosen_winv','winchosen_lossv',...
    'lossnotchosen_bothv','lossnotchosen_winv','lossnotchosen_lossv',...
    'winvslosschosen_bothv','winvslosschosen_winv','winvslosschosen_lossv',});
T_alphas=array2table([switchprob_win.session1.all.all,switchprob_nowin.session1.all.all,...
    switchprob_loss.session1.all.all,switchprob_noloss.session1.all.all],...
    'VariableNames',{'switchprob_win_bothv','switchprob_win_winv','switchprob_win_lossv',...
                                  'switchprob_nowin_bothv','switchprob_nowin_winv','switchprob_nowin_lossv',...
                                  'switchprob_loss_bothv','switchprob_loss_winv','switchprob_loss_lossv',...
                                  'switchprob_noloss_bothv','switchprob_noloss_winv','switchprob_noloss_lossv'});
T_ques=array2table([QIDS.session1.all.all,sSTAI.session1.all.all,tSTAI.session1.all.all,SHAPS.session1.all.all,...
    DAQ.session1.all.all,RRS.session1.all.all],...
    'VariableNames',{'QIDS','sSTAI','tSTAI','SHAPS','DAQ','RRS'});
T_2lr_2b=array2table([alpha_rew_2lr_2b.session1.all.all,alpha_loss_2lr_2b.session1.all.all,...
    beta_rew_2lr_2b.session1.all.all,beta_loss_2lr_2b.session1.all.all],...
    'VariableNames',{'rew_alpha_2lr_2b_bothv','rew_alpha_2lr_2b_winv','rew_alpha_2lr_2b_lossv',...
                                  'loss_alpha_2lr_2b_bothv','loss_alpha_2lr_2b_winv','loss_alpha_2lr_2b_lossv',...
                                  'rew_beta_2lr_2b_bothv','rew_beta_2lr_2b_winv','rew_beta_2lr_2b_lossv',...
                                  'loss_beta_2lr_2b_bothv','loss_beta_2lr_2b_winv','loss_beta_2lr_2b_lossv'});
T_PNPE_2lr_2b=array2table([alpha_PPE_2lr_2b.session1.all.all,alpha_NPE_2lr_2b.session1.all.all,...
    beta_rew_PNPE_2lr_2b.session1.all.all,beta_loss_PNPE_2lr_2b.session1.all.all],...
    'VariableNames',{'Positive_alpha_PNPE_2lr_2b_bothv','Positive_alpha_PNPE_2lr_2b_winv','Positive_alpha_PNPE_2lr_2b_lossv',...
                                  'Negative_alpha_PNPE_2lr_2b_bothv','Negative_alpha_PNPE_2lr_2b_winv','Negative_alpha_PNPE_2lr_2b_lossv',...
                                  'rew_beta_PNPE_2lr_2b_bothv','rew_beta_PNPE_2lr_2b_winv','rew_beta_PNPE_2lr_2b_lossv',...
                                  'loss_beta_PNPE_2lr_2b_bothv','loss_beta_PNPE_2lr_2b_winv','loss_beta_PNPE_2lr_2b_lossv'});
T_2lr_1b=array2table([alpha_rew_2lr_1b.session1.all.all,alpha_loss_2lr_1b.session1.all.all,beta_2lr_1b.session1.all.all],...
    'VariableNames',{'rew_alpha_2lr_1b_bothv','rew_alpha_2lr_1b_winv','rew_alpha_2lr_1b_lossv',...
                                  'loss_alpha_2lr_1b_bothv','loss_alpha_2lr_1b_winv','loss_alpha_2lr_1b_lossv',...
                                  'beta_2lr_1b_bothv','beta_2lr_1b_winv','beta_2lr_1b_lossv'});
T_2lr_1b_bias=array2table([alpha_rew_2lr_1b_bias.session1.all.all,alpha_loss_2lr_1b_bias.session1.all.all,beta_2lr_1b_bias.session1.all.all],...
    'VariableNames',{'rew_alpha_2lr_1b_bias_bothv','rew_alpha_2lr_1b_bias_winv','rew_alpha_2lr_1b_bias_lossv',...
                                  'loss_alpha_2lr_1b_bias_bothv','loss_alpha_2lr_1b_bias_winv','loss_alpha_2lr_1b_bias_lossv',...
                                  'beta_2lr_1b_bias_bothv','beta_2lr_1b_bias_winv','beta_2lr_1b_bias_lossv'});
T_1lr_1b=array2table([alpha_1lr_1b.session1.all.all,beta_1lr_1b.session1.all.all],...
    'VariableNames',{'alpha_1lr_1b_bothv','alpha_1lr_1b_winv','alpha_1lr_1b_lossv',...
                                  'beta_1lr_1b_bothv','beta_1lr_1b_winv','beta_1lr_1b_lossv'});
T=[T_s,T_group,T_treatment,T_ques,T_performance,T_alphas,T_2lr_2b,T_2lr_1b,T_2lr_1b_bias,T_PNPE_2lr_2b,T_1lr_1b];
writetable(T,[datadir,'2017_clinical_study_session1'])

% for group={'reboxetine','placebo','all'}
%     for v=1:2
% plot_bargarph_BICs
%     end
% end
% %% AIC
% for i=1:size(result_2lr_2b,1)
%     for j=1:4%size(result_2lr_2b,2)
%           AIC_2lr_2b(v,i,j)=getfield(result_2lr_2b,{i,j+1},'AIC');
%           AIC_2lr_1b(v,i,j)=getfield(result_2lr_1b,{i,j+1},'AIC');
%           AIC_1lr_1b(v,i,j)=getfield(result_1lr_1b,{i,j+1},'AIC');
%           AIC_1lr_1b_opt1(v,i,j)=getfield(result_PNPE_2lr_1b_opt1,{i,j+1},'AIC');
%           AIC_PNPE_2lr_2b(v,i,j)=getfield(result_PNPE_2lr_2b,{i,j+1},'AIC');
%           AIC_PNPE_2lr_1b(v,i,j)=getfield(result_PNPE_2lr_1b,{i,j+1},'AIC');
%           AIC_PNPE_2lr_1b_opt1(v,i,j)=getfield(result_PNPE_2lr_1b_opt1,{i,j+1},'AIC');
%           AIC_2lr_1b_sdata(i).choice_stickness(v,i,j)=getfield(result_2lr_1b_sdata(i).choice_stickness,{i,j+1},'AIC');
%     end
% end
% plot_bargarph_AICs
% 
% %%
% T_s=table(Ques.subnum,Ques.treatment,'VariableNames',{'subjectnum','treatment'});
% T_rew_alphas=array2table([squeeze(inv_logit(alpha_rew_2lr_1b.all(1,:,:))),squeeze(inv_logit(alpha_rew_2lr_1b.all(2,:,:)))],...
%     'VariableNames',{'rew_alpha_2lr1b_visit1_bothv','rew_alpha_2lr1b_visit1_winv','rew_alpha_2lr1b_visit1_lossv','rew_alpha_2lr1b_visit1_boths',...
%                                     'rew_alpha_2lr1b_visit2_bothv','rew_alpha_2lr1b_visit2_winv','rew_alpha_2lr1b_visit2_lossv','rew_alpha_2lr1b_visit2_boths'});
% T_loss_alphas=array2table([squeeze(inv_logit(alpha_loss_2lr_1b.all(1,:,:))),squeeze(inv_logit(alpha_loss_2lr_1b.all(2,:,:)))],...
%     'VariableNames',{'loss_alpha_2lr1b_visit1_bothv','loss_alpha_2lr1b_visit1_winv','loss_alpha_2lr1b_visit1_lossv','loss_alpha_2lr1b_visit1_boths',...
%                                     'loss_alpha_2lr1b_visit2_bothv','loss_alpha_2lr1b_visit2_winv','loss_alpha_2lr1b_visit2_lossv','loss_alpha_2lr1b_visit2_boths'});
% T_betas=array2table([squeeze(log(beta_2lr_1b.all(1,:,:))),squeeze(log(beta_2lr_1b.all(2,:,:)))],...
%     'VariableNames',{'beta_2lr1b_visit1_bothv','beta_2lr1b_visit1_winv','beta_2lr1b_visit1_lossv','beta_2lr1b_visit1_boths',...
%                                     'beta_2lr1b_visit2_bothv','beta_2lr1b_visit2_winv','beta_2lr1b_visit2_lossv','beta_2lr1b_visit2_boths'});  
% %rm anova model for alphas
% T_stat_alphas=[T_s,T_rew_alphas,T_loss_alphas];
% within=table(repelem({'win';'loss'},8,1),repmat(repelem({'visit1';'visit2'},4,1),2,1),repmat({'winv';'winv';'wins';'wins'},4,1),repmat({'lossv';'losss';'lossv';'losss'},4,1),...
%      'VariableNames',{'valence','visit','winvol','lossvol'});
%  rm=fitrm(T_stat_alphas,'rew_alpha_2lr1b_visit1_bothv-loss_alpha_2lr1b_visit2_boths ~ treatment','WithinDesign',within);
%  ranovatbl_alpha_2lr_1b=ranova(rm,'WithinModel','valence*visit*winvol*lossvol')
% %rm anova model for batas
% T_stat_betas=[T_s,T_betas];
% within=table(repelem({'visit1';'visit2'},4,1),repmat({'winv';'winv';'wins';'wins'},2,1),repmat({'lossv';'losss';'lossv';'losss'},2,1),...
%      'VariableNames',{'visit','winvol','lossvol'});
% rm=fitrm(T_stat_betas,'beta_2lr1b_visit1_bothv-beta_2lr1b_visit2_boths ~ treatment','WithinDesign',within);
% ranovatbl_beta_2lr1b=ranova(rm,'WithinModel','visit*winvol*lossvol')
% T=[T_s, T_rew_alphas,T_loss_alphas,T_betas];
%  writetable(T,[datadir,'2lr1b_model'])
% %% do stats for 2lr1b bias model
% T_s=table(Ques.subnum,Ques.treatment,'VariableNames',{'subjectnum','treatment'});
% T_rew_alphas=array2table([squeeze(inv_logit(alpha_rew_2lr_1b_bias.all(1,:,:))),squeeze(inv_logit(alpha_rew_2lr_1b_bias.all(2,:,:)))],...
%     'VariableNames',{'rew_alpha_2lr1b_bias_visit1_bothv','rew_alpha_2lr1b_bias_visit1_winv','rew_alpha_2lr1b_bias_visit1_lossv','rew_alpha_2lr1b_bias_visit1_boths',...
%                                     'rew_alpha_2lr1b_bias_visit2_bothv','rew_alpha_2lr1b_bias_visit2_winv','rew_alpha_2lr1b_bias_visit2_lossv','rew_alpha_2lr1b_bias_visit2_boths'});
% T_loss_alphas=array2table([squeeze(inv_logit(alpha_loss_2lr_1b_bias.all(1,:,:))),squeeze(inv_logit(alpha_loss_2lr_1b_bias.all(2,:,:)))],...
%     'VariableNames',{'loss_alpha_2lr1b_bias_visit1_bothv','loss_alpha_2lr1b_bias_visit1_winv','loss_alpha_2lr1b_bias_visit1_lossv','loss_alpha_2lr1b_bias_visit1_boths',...
%                                     'loss_alpha_2lr1b_bias_visit2_bothv','loss_alpha_2lr1b_bias_visit2_winv','loss_alpha_2lr1b_bias_visit2_lossv','loss_alpha_2lr1b_bias_visit2_boths'});
% T_betas=array2table([squeeze(log(beta_2lr_1b_bias.all(1,:,:))),squeeze(log(beta_2lr_1b_bias.all(2,:,:)))],...
%     'VariableNames',{'beta_2lr1b_bias_visit1_bothv','beta_2lr1b_bias_visit1_winv','beta_2lr1b_bias_visit1_lossv','beta_2lr1b_bias_visit1_boths',...
%                                     'beta_2lr1b_bias_visit2_bothv','beta_2lr1b_bias_visit2_winv','beta_2lr1b_bias_visit2_lossv','beta_2lr1b_bias_visit2_boths'});
% T_biases=array2table([squeeze(abs(bias_2lr_1b_bias.all(1,:,:))),squeeze(abs(bias_2lr_1b_bias.all(2,:,:)))],...
%     'VariableNames',{'bias_2lr1b_bias_visit1_bothv','bias_2lr1b_bias_visit1_winv','bias_2lr1b_bias_visit1_lossv','bias_2lr1b_bias_visit1_boths',...
%                                     'bias_2lr1b_bias_visit2_bothv','bias_2lr1b_bias_visit2_winv','bias_2lr1b_bias_visit2_lossv','bias_2lr1b_bias_visit2_boths'}); 
%                                 
% %rm anova model for alphas
% T_stat_alphas=[T_s,T_rew_alphas,T_loss_alphas];
% within=table(repelem({'win';'loss'},8,1),repmat(repelem({'visit1';'visit2'},4,1),2,1),repmat({'winv';'winv';'wins';'wins'},4,1),repmat({'lossv';'losss';'lossv';'losss'},4,1),...
%      'VariableNames',{'valence','visit','winvol','lossvol'});
%  rm=fitrm(T_stat_alphas,'rew_alpha_2lr1b_bias_visit1_bothv-loss_alpha_2lr1b_bias_visit2_boths ~ treatment','WithinDesign',within);
%  ranovatbl_alphas_2lr1b_bias=ranova(rm,'WithinModel','valence*visit*winvol*lossvol')
% %rm anova model for batas
% T_stat_betas=[T_s,T_betas];
% within=table(repelem({'visit1';'visit2'},4,1),repmat({'winv';'winv';'wins';'wins'},2,1),repmat({'lossv';'losss';'lossv';'losss'},2,1),...
%      'VariableNames',{'visit','winvol','lossvol'});
% rm=fitrm(T_stat_betas,'beta_2lr1b_bias_visit1_bothv-beta_2lr1b_bias_visit2_boths ~ treatment','WithinDesign',within);
% ranovatbl_betas_2lr1b_bias=ranova(rm,'WithinModel','visit*winvol*lossvol')
% %rm anova model for bias term
% T_stat_biases=[T_s,T_biases];
% within=table(repelem({'visit1';'visit2'},4,1),repmat({'winv';'winv';'wins';'wins'},2,1),repmat({'lossv';'losss';'lossv';'losss'},2,1),...
%      'VariableNames',{'visit','winvol','lossvol'});
%  rm=fitrm(T_stat_biases,'bias_2lr1b_bias_visit1_bothv-bias_2lr1b_bias_visit2_boths ~ treatment','WithinDesign',within);
%  ranovatbl_bias_2lr1b_bias=ranova(rm,'WithinModel','visit*winvol*lossvol')
%  
%  % write table to txt for spss
%  T=[T_s, T_rew_alphas,T_loss_alphas,T_betas,T_biases];
%  writetable(T,[datadir,'2lr1b_bias_model'])

%% back up
for session={'session1','session2','session3'}
    for ss=1:size(subspool,2)
    oridata=importdata([datadir,num2str(subspool(ss)),'/',num2str(subspool(ss)),'_',strrep(session{1},'n','n_'),'_vol_train_log.dat']);
    sort_data_into_blocks
        %fit models
        for i=1:blkn
        sdata(i).choice(sdata(i).choice==2)=0;
        %2alpha2beta
       % result_2lr_2b.(session{1}).all.all(ss,i,:)=fit_linked_2lr_2beta_il_rev(sdata(i).information,sdata(i).choice, start, abandontn);
        %1alpha1beta
        %result_1lr_1b.(session{1}).all.all(ss,i,:)=fit_linked_1lr_1beta_il(sdata(i).information,sdata(i).choice, start, abandontn);
        %2alpha1beta
        %result_2lr_1b.(session{1}).all.all(ss,i,:)=fit_linked_2lr_1beta_il(sdata(i).information,sdata(i).choice, start, abandontn);
        %2alphas1beta+1bias term
        result_2lr_1b_bias.(session{1}).all.all(ss,i,:)=fit_linked_2lr_beta_add(sdata(i).information,sdata(i).choice, start, abandontn);
        %2alpha1beta+sdata(i).choice_stickness
      %  result_2lr_1b_choice_stickness.(session{1}).all.all(ss,i,:)=fit_linked_2lr_1beta_il_choice_stickness(sdata(i).information,sdata(i).choice, start, abandontn);
        %opt1result positive negative PE learning rates
      %  result_PNPE_2lr_1b_opt1.(session{1}).all.all(ss,i,:)=fit_linked_PNPE_2lr_1beta_opt1results(sdata(i).information,sdata(i).choice, 0, abandontn);
              end
    end
end
%%
clear all

load('/home/wlin/Documents/2017_clinical_study/code/pupil/elife_analysis/pupil.mat')
comsubs=ismember(sublist.all,session1_sublist.all);
newgroup=[ones(size(session1_sublist.G1,2),1);ones(size(session1_sublist.G2,2),1)*2;ones(size(session1_sublist.G3,2),1)*3];

[rwin,pwin]=corr(mean(adj_means_rew_undelivered.session_1.all.all(:,1001:end),2), delta_win.session1.all.all(comsubs));
anova1(mean(adj_pupil_delta_win.session_1.all.all(:,1001:1500),2),newgroup)
x=3501;
anova1(mean(adj_pupil_delta_win.session_1.all.all(:,x:x+500),2),newgroup)
[rloss,ploss]=corr(mean(adj_means_pun_delivered.session_1.all.all(:,1001:end),2), delta_loss.session1.all.all(comsubs));
[rloss,ploss]=corr(mean(adj_pupil_delta_loss.session_1.all.all(:,1001:end),2), delta_loss.session1.all.all(comsubs));
scatterplot_3groups(mean(adj_pupil_delta_loss.session_1.all.all(:,1001:end),2),delta_loss.session1.all.all(comsubs),...
    '-','-','tmp','Pearson',newgroup,figdir)
comsubs=ismember(sublist.G3,session1_sublist.G3);
[rwin,pwin]=corr(mean(adj_means_pun_delivered.session_1.G3.all(:,1001:end),2), delta_loss.session1.G3.all(comsubs))

[rwin,pwin]=corr(mean(means_rew_undelivered.session_1.all.all(:,1001:end),2), delta_win.session1.all.all(comsubs));
%%
blk=3;
[rloss,ploss]=corr(squeeze(means_pun_delivered.session_1.all.all(:,blk,:)),switchprob_loss.session1.all.all(comsubs,blk));
plot(rloss)
[rloss,ploss]=corr(squeeze(mean(means_pun_delivered.session_1.all.all(:,blk,1001:2000),3)),switchprob_loss.session1.all.all(comsubs,blk))
%%
blk=3;%also significant when blk=1;
[rloss,ploss]=corr(squeeze(means_pun_undelivered.session_1.all.all(:,blk,:)),switchprob_noloss.session1.all.all(comsubs,blk));
plot(rloss)
[rloss,ploss]=corr(squeeze(mean(means_pun_undelivered.session_1.all.all(:,blk,1001:2000),3)),switchprob_noloss.session1.all.all(comsubs,blk))

%%
blk=2;
[rwin,pwin]=corr(squeeze(means_rew_delivered.session_1.all.all(:,blk,:)),switchprob_win.session1.all.all(comsubs,blk));
plot(rwin)
[rwin,pwin]=corr(squeeze(mean(means_rew_delivered.session_1.all.all(:,blk,1001:2000),3)),switchprob_win.session1.all.all(comsubs,blk))
%%
blk=3;%???why
[rwin,pwin]=corr(squeeze(means_rew_undelivered.session_1.all.all(:,blk,:)),switchprob_nowin.session1.all.all(comsubs,blk));
plot(rwin)
[rwin,pwin]=corr(squeeze(mean(means_rew_undelivered.session_1.all.all(:,blk,1001:2000),3)),switchprob_nowin.session1.all.all(comsubs,blk))

%%
[rwin,pwin]=corr(adj_means_pun_delivered.session_1.all.all,switchprob_loss.session1.all.all(comsubs,3)-switchprob_loss.session1.all.all(comsubs,2));
plot(rwin)
[rwin,pwin]=corr(adj_means_pun_undelivered.session_1.all.all,switchprob_noloss.session1.all.all(comsubs,3)-switchprob_noloss.session1.all.all(comsubs,2));
plot(rwin)