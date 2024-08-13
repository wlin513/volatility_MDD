%% analyse pupil timecourse model free based on outcome delivered
%% 
clear all;
clc;
%%
cd('/home/wlin/Documents/2017_clinical_study/code/pupil/elife_analysis');
getfolders;
addpath('/home/wlin/Documents/2017_clinical_study/code/plot_code/')
addpath('/home/wlin/Documents/2017_clinical_study/code/useful_fun/')
%% load subject data pool
session1_sublist.G1=1001:1041;session1_sublist.G1(session1_sublist.G1==1036)=[];
session1_sublist.G2=2001:2042; session1_sublist.G2(session1_sublist.G2==2033)=[];session1_sublist.G2(session1_sublist.G2==2019)=[];
session1_sublist.G3=3001:3041; session1_sublist.G3(session1_sublist.G3==3033)=[];
session1_sublist.G1(session1_sublist.G1==1004)=[];%exclueded from beh analysis
session1_sublist.G1(session1_sublist.G1==1039)=[];%exclueded from beh analysis
session1_sublist.G2(session1_sublist.G1==2034)=[];%exclueded from beh analysis
% exclude participants for session 1
session1_sublist.G1(session1_sublist.G1==1005)=[];%no session 1 pupil data
%interpolated more than 60% of the pupil data
session1_sublist.G1(session1_sublist.G1==1031)=[];% s1 62% s3 72%
session1_sublist.G3(session1_sublist.G3==3024)=[];%somehow all session1 results are nans
session1_sublist.G3(session1_sublist.G3==3009)=[];%s1 65% s3 65% 
session1_sublist.G3(session1_sublist.G3==3029)=[];%s1 72% s3 75% 
session1_sublist.G3(sublist.G3==3040)=[];%couldn't read event asc file but edf2mat file is fine. double check when have time/needed

% % exclusion: exclude participant who some of the data were missing
% %G1
% sublist.G1(sublist.G1==1009)=[];%no session 2 pupil data
% sublist.G1(sublist.G1==1014)=[];%no session 3 pupil data
% %G2
% sublist.G2(sublist.G2==2012)=[];%no session 2 pupil data
% sublist.G2(sublist.G2==2016)=[];%no session 2 pupil data
% %G3
% sublist.G1(sublist.G1==1016)=[];% s3 71%
% sublist.G2(sublist.G2==2025)=[];% s2 80% s3 61%
% sublist.G2(sublist.G2==2008)=[];% s3 85%
% sublist.G2(sublist.G2==2024)=[];% s3 88%
% sublist.G2(sublist.G2==2039)=[];% s3 64%
% sublist.G2(sublist.G2==2038)=[];% s2 75%
% sublist.G3(sublist.G3==3006)=[];%s3 66%
% sublist.G3(sublist.G3==3025)=[];%s2 65%
% sublist.G3(sublist.G3==3028)=[];%s2 88% s3 70%
% sublist.G3(sublist.G3==3041)=[];%s2 75%


session1_sublist.all=[session1_sublist.G1,session1_sublist.G2,session1_sublist.G3];
group=[ones(size(session1_sublist.G1,2),1);ones(size(session1_sublist.G2,2),1)*2;ones(size(session1_sublist.G3,2),1)*3];
%% get training info
load([datadir,'finaldata.mat'])
for subset={'G1','G2','G3','all'}
   subspool=session1_sublist.(subset{1});
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
QIDS=get_p1vital_info(dataout,session1_sublist,trainingtype,'QIDS');
sSTAI=get_p1vital_info(dataout,session1_sublist,trainingtype,'sSTAI');
tSTAI=get_p1vital_info(dataout,session1_sublist,trainingtype,'tSTAI');
DAQ=get_p1vital_info(dataout,session1_sublist,trainingtype,'DAQ');
RRS=get_p1vital_info(dataout,session1_sublist,trainingtype,'RRS');

%% experiment information setting
tn=80;
nblk=3;

blkname={'both volatile','win volatile','loss volatile'};
blockorders=[1,3,2;1,2,3];%1 for both volatile 2 for win volatile 3 for loss volatile
abandontn=10;%abandon first 10 trails for each block
%%  

    ssublist=session1_sublist.all;
    for session={'session_1'}%,'session_2','session_3'}
       for sub=1:size(ssublist,2)    
                name=num2str(ssublist(sub));
                sdatadir=[datadir,name,'/'];     
                %% load pupil data 
                file_name=['out_',name,'_',session{1},'_normalised_amongst_3sessions.mat'];
                load([sdatadir,file_name])
                load([sdatadir,name,'_',session{1},'_pupil_exclude_trials.mat']);
                rew_include=1-rew_exclude'; % to exclude trials where more than 50% of the data were missing
                pun_include=1-pun_exclude';
                rew_data=out.rew_outcome;
                pun_data=out.pun_outcome;
                %% loads the behavior data
                oridata=importdata([sdatadir,name,'_',session{1},'_vol_train_log.dat']);
                data=struct;
                for i=1:length(oridata.colheaders)
                    data.(oridata.colheaders{i})=oridata.data(:,i);
                end   
                data.version=str2double(extractAfter(oridata.textdata{1},'order_num'));
                ifchosen_rew=data.shape1_win;
                ifchosen_rew(data.button_pressed==2)=data.shape2_win(data.button_pressed==2);
                ifchosen_pun=data.shape1_loss;
                ifchosen_pun(data.button_pressed==2)=data.shape2_loss(data.button_pressed==2);

                rew_order=1-data.outcome_order;%1 for win first 0 for loss first
                pun_order=data.outcome_order;%1 for loss first 0 for win first
                %% baseline corrections 
                rewbaselines=mean(rew_data(:,501:1000),2); %every 500 represents 1 seconds
                rew_data=rew_data-rewbaselines;
                punbaselines=mean(pun_data(:,501:1000),2); 
                pun_data=pun_data-punbaselines;           
                %% exclude trails and then calculate mean across trials for each trial type
                includetrials=repmat([false(abandontn,1);true(tn-abandontn,1)],nblk,1);
                rewnantrials=~isnan(mean(rew_data,2));
                punnantrials=~isnan(mean(pun_data,2));
                summary_nantrials.(session{1})(sub,:)=[sum(isnan(mean(rew_data,2))),sum(isnan(mean(pun_data,2)))];
                % generate block index so that every subject's data is organised in the order of 1 both vol,2 win vol, 3 loss vol
                blkorder=blockorders(data.version,:);
                for blk=1:nblk   
                        blkindex=false(size(rew_data,1),1);
                        blkloc=find(blkorder==blk);
                        blkindex((blkloc-1)*tn+1:blkloc*tn)=true;
                        rew_delivered=rew_data(ifchosen_rew==1&blkindex&includetrials&rewnantrials&rew_include,:);
                        rew_undelivered=rew_data(ifchosen_rew==0&blkindex&includetrials&rewnantrials&rew_include,:);
                        means_rew_delivered.(session{1}).all.all(sub,blk,:)=mean(rew_delivered,1); %take means for each timepoint
                        means_rew_undelivered.(session{1}).all.all(sub,blk,:)=mean(rew_undelivered,1);
                        pun_delivered=pun_data(ifchosen_pun==1&blkindex&includetrials&punnantrials&pun_include,:);
                        pun_undelivered=pun_data(ifchosen_pun==0&blkindex&includetrials&punnantrials&pun_include,:);
                        means_pun_delivered.(session{1}).all.all(sub,blk,:)=mean(pun_delivered,1);
                        means_pun_undelivered.(session{1}).all.all(sub,blk,:)=mean(pun_undelivered,1);                       

                        rew_bslines=rewbaselines(blkindex&includetrials&rewnantrials&rew_include);
                        means_rew_bslines.(session{1}).all.all(sub,blk)=mean(rew_bslines,1);
                        pun_bslines=punbaselines(blkindex&includetrials&punnantrials&pun_include);
                        means_pun_bslines.(session{1}).all.all(sub,blk)=mean(pun_bslines,1);
                end
                %%mean for each conditions
                %mean for win vol blocks
                blkindex=false(size(rew_data,1),1);
                blkloc=find(blkorder==1|blkorder==2);
                for i=1:length(blkloc)
                    blkindex((blkloc(i)-1)*tn+1:blkloc(i)*tn)=true;
                end
                rew_delivered=rew_data(ifchosen_rew==1&blkindex&includetrials&rewnantrials&rew_include,:);
                rew_undelivered=rew_data(ifchosen_rew==0&blkindex&includetrials&rewnantrials&rew_include,:);
                delta_win_vol.(session{1}).all.all(sub,:)=mean(rew_delivered,1)-mean(rew_undelivered,1);
                
               %mean for loss vol blocks
                blkindex=false(size(rew_data,1),1);
                blkloc=find(blkorder==1|blkorder==3);
                for i=1:length(blkloc)
                    blkindex((blkloc(i)-1)*tn+1:blkloc(i)*tn)=true;
                end
                pun_delivered=pun_data(ifchosen_pun==1&blkindex&includetrials&punnantrials&pun_include,:);
                pun_undelivered=pun_data(ifchosen_pun==0&blkindex&includetrials&punnantrials&pun_include,:);
                delta_loss_vol.(session{1}).all.all(sub,:)=mean(pun_delivered,1)-mean(pun_undelivered,1);
             
        end
    
    delta_win.(session{1}).all.all=means_rew_delivered.(session{1}).all.all-means_rew_undelivered.(session{1}).all.all; %pupil response to rewards when chosen vs unchosen rewarded
    delta_loss.(session{1}).all.all=means_pun_delivered.(session{1}).all.all-means_pun_undelivered.(session{1}).all.all;
    delta_winvsloss.(session{1}).all.all=delta_win.(session{1}).all.all-delta_loss.(session{1}).all.all;
    adj_delta_win.(session{1}).all.all=squeeze(delta_win.(session{1}).all.all(:,2,:)-delta_win.(session{1}).all.all(:,3,:));
    adj_delta_loss.(session{1}).all.all=squeeze(delta_loss.(session{1}).all.all(:,3,:)-delta_loss.(session{1}).all.all(:,2,:));
    adj_means_rew_delivered.(session{1}).all.all=squeeze(means_rew_delivered.(session{1}).all.all(:,2,:)-means_rew_delivered.(session{1}).all.all(:,3,:));
    adj_means_pun_delivered.(session{1}).all.all=squeeze(means_pun_delivered.(session{1}).all.all(:,3,:)-means_pun_delivered.(session{1}).all.all(:,2,:));
    adj_means_rew_undelivered.(session{1}).all.all=squeeze(means_rew_undelivered.(session{1}).all.all(:,2,:)-means_rew_undelivered.(session{1}).all.all(:,3,:));
    adj_means_pun_undelivered.(session{1}).all.all=squeeze(means_pun_undelivered.(session{1}).all.all(:,3,:)-means_pun_undelivered.(session{1}).all.all(:,2,:));
    end
    delta_win=get_full_structure(delta_win,trainingtype,session1_sublist);
    delta_loss=get_full_structure(delta_loss,trainingtype,session1_sublist);
    delta_winvsloss=get_full_structure(delta_winvsloss,trainingtype,session1_sublist);
    adj_delta_win=get_full_structure(adj_delta_win,trainingtype,session1_sublist);
    adj_delta_loss=get_full_structure(adj_delta_loss,trainingtype,session1_sublist);
    means_rew_delivered=get_full_structure(means_rew_delivered,trainingtype,session1_sublist);
    means_rew_undelivered=get_full_structure(means_rew_undelivered,trainingtype,session1_sublist);
    means_pun_delivered=get_full_structure(means_pun_delivered,trainingtype,session1_sublist);
    means_pun_undelivered=get_full_structure(means_pun_undelivered,trainingtype,session1_sublist);
    adj_means_rew_delivered=get_full_structure(adj_means_rew_delivered,trainingtype,session1_sublist);
    adj_means_rew_undelivered=get_full_structure(adj_means_rew_undelivered,trainingtype,session1_sublist);
    adj_means_pun_delivered=get_full_structure(adj_means_pun_delivered,trainingtype,session1_sublist);
    adj_means_pun_undelivered=get_full_structure(adj_means_pun_undelivered,trainingtype,session1_sublist);
%% plot figures
save('pupil.mat')
plot_pupil_for_each_block(delta_win,delta_loss,'session_1','received vs not',blkname,figdir)
plot_pupil_for_each_block(delta_winvsloss,delta_winvsloss,'session_1','win vs loss',blkname,figdir)
plot_pupil_for_each_block(means_rew_delivered,means_pun_delivered,'session_1','received')
plot_pupil_for_each_block(means_rew_undelivered,means_pun_undelivered,'session_1','not received')
plot_pupil_for_each_block(delta_win,delta_loss,'session_2','received vs not',blkname,figdir)
plot_pupil_for_mean(adj_delta_win,adj_delta_loss,'session_1','vol vs stable',figdir)

plot_pupil_for_mean(adj_means_rew_delivered,adj_means_pun_delivered,'session_1','vol vs stable',figdir)
plot_pupil_for_mean(adj_means_rew_undelivered,adj_means_pun_undelivered,'session_1','vol vs stable',figdir)

% %%
% plot_pupil_vol_vs_stable
% plot_pupil_received_each_block
% plot_pupil_not_received_each_block
% plot_pupil_each_block
% barplot_baselines
% plot_pupil_received_and_not_for_winv_and_lossv_block
% plot_pupil_received_and_not_for_semi_vs_balanced_blocks
% plot_pupil_response_diff_vol_sta_to_otherschedule