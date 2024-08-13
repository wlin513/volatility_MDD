%% 
clear all;
clc;
close all;
%%
cd('D:/2017_clinical_study/code');
getfolders;
addpath('plot_code/')
addpath('useful_fun/')
%% load subject data pool

sublist.G1=1001:1041;
sublist.G2=2001:2042; sublist.G2(sublist.G2==2033)=[];
sublist.G3=3001:3041;
% exclusion: exclude participant who chose shape 1 or left side ???
%G1
sublist.G1(sublist.G1==1004)=[];%who chose shape 1 for 66% 0% 97% for the 3 blocks respectively in session1
sublist.G1(sublist.G1==1039)=[];%who chose left side for 38.75% 0% 42.5% for the 3 blocks respectively in session1
% sublist.G1(sublist.G1==1017)=[];%who chose shape 1 for 45% 87.5% 1.25% for the 3 blocks respectively in session 1
%sublist.G1(sublist.G1==1010)=[];%who chose shape 1 for 32.5% 85% 95% for the 3 blocks respectively in session 1
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
DAQ=get_p1vital_info(dataout,sublist,trainingtype,'DAQ');
RRS=get_p1vital_info(dataout,sublist,trainingtype,'RRS');
SHAPS=get_p1vital_info(dataout,sublist,trainingtype,'SHAPS');

%% experiment information setting
tn=80;
blkn=3;

blkname={'both volatile','win volatile','loss volatile'};
blkindex=[1,3,2;1,2,3];%1 for both volatile 2 for win volatile 3 for loss volatile
resp_made=ones(tn,1);
%%  calculate switch probabilities
for session={'session1'}%,'session2','session3'}
subspool=sublist.all;
%
for ss=1:size(subspool,2)
%
oridata=importdata([datadir,num2str(subspool(ss)),'/',num2str(subspool(ss)),'_',strrep(session{1},'n','n_'),'_vol_train_log.dat']);

sort_data_into_blocks

%checking probability of chossing left or shape A in each block
for i=1:blkn
    tmp_pct_leftchoice(i,:)=mean(sdata(i).choiceside);
    tmp_pct_shape1choice(i,:)=mean(sdata(i).choice-1);
    
    %percentages of winchosen and loss not chosen
    tmp_pct_winchosen(i,:)=mean(sdata(i).winchosen);
    tmp_pct_losschosen(i,:)=mean(sdata(i).losschosen);

    %calculate swtich probabilty for each block
    tmp_switchprob(i,:)=cal_switch_prob(sdata(i).information,sdata(i).choice,resp_made);
end
pct_leftchoice.(session{1}).all.all(ss,:)=tmp_pct_leftchoice;
pct_shape1choice.(session{1}).all.all(ss,:)=tmp_pct_shape1choice;
pct_winchosen.(session{1}).all.all(ss,:)=tmp_pct_winchosen;
pct_lossnotchosen.(session{1}).all.all(ss,:)=1-tmp_pct_losschosen;
tmp_pct_leftchoice=[];
tmp_pct_shape1choice=[];
tmp_pct_winchosen=[];
tmp_pct_losschosen=[];
switchprob.(session{1}).all.all(ss,:)=tmp_switchprob;
clear tmp_switchprob
end
%percentages of win vs loss chosen
pct_winvsloss_chosen.(session{1}).all.all=pct_winchosen.(session{1}).all.all-(1-pct_lossnotchosen.(session{1}).all.all);
end

pct_leftchoice=get_full_structure(pct_leftchoice,trainingtype,sublist);
pct_winvsloss=get_full_structure(pct_winvsloss_chosen,trainingtype,sublist);
pct_shape1choice=get_full_structure(pct_shape1choice,trainingtype,sublist);
switchprob=get_full_structure(switchprob,trainingtype,sublist);
pct_winchosen=get_full_structure(pct_winchosen,trainingtype,sublist);
pct_lossnotchosen=get_full_structure(pct_lossnotchosen,trainingtype,sublist);
%%
save('model_free_results_session1.mat')
%% repeated anova
T_s=array2table(sublist.all','VariableNames',{'subject'});
T_group=array2table(groups,'VariableNames',{'group'});
T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});
T_performance=array2table([pct_winchosen.session1.all.all,pct_lossnotchosen.session1.all.all,pct_winvsloss.session1.all.all],...
    'VariableNames',{'winchosen_bothv','winchosen_winv','winchosen_lossv',...
    'lossnotchosen_bothv','lossnotchosen_winv','lossnotchosen_lossv',...
    'winvslosschosen_bothv','winvslosschosen_winv','winvslosschosen_lossv',});
T=[T_s,T_group,T_treatment,T_performance];
% within=table(repelem({'win';'loss'},3,1),repmat({'winv';'winv';'wins'},2,1),repmat({'lossv';'losss';'lossv'},2,1),...
%      'VariableNames',{'valence','winvol','lossvol'});
within=table(repelem({'win';'loss'},3,1),repmat({'1';'2';'3'},2,1),...
     'VariableNames',{'valence','blktype'});
rm=fitrm(T,'winchosen_bothv-lossnotchosen_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','valence*blktype');

%
T_performance=array2table([pct_winchosen.session1.all.all(:,2:3),pct_lossnotchosen.session1.all.all(:,2:3),pct_winvsloss.session1.all.all(:,2:3)],...
    'VariableNames',{'winchosen_winv','winchosen_lossv',...
    'lossnotchosen_winv','lossnotchosen_lossv',...
    'winvslosschosen_winv','winvslosschosen_lossv',});
T=[T_s,T_group,T_treatment,T_performance];
% within=table(repelem({'win';'loss'},3,1),repmat({'winv';'winv';'wins'},2,1),repmat({'lossv';'losss';'lossv'},2,1),...
%      'VariableNames',{'valence','winvol','lossvol'});
within=table(repelem({'win';'loss'},2,1),repmat({'winv';'lossv'},2,1),...
     'VariableNames',{'valence','blktype'});
rm=fitrm(T,'winchosen_winv-lossnotchosen_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','valence*blktype');
%
within=table({'1';'2';'3'},...
     'VariableNames',{'blktype'});
rm=fitrm(T,'winchosen_bothv-winchosen_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','blktype');
%
within=table({'1';'2';'3'},...
     'VariableNames',{'blktype'});
rm=fitrm(T,'lossnotchosen_bothv-lossnotchosen_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','blktype');
%plot the effect
plot_bargarph_3group_comparison(pct_winchosen,'pct-winchosen','session1','percentage chosen',[0.4,0.7],blkname,figdir)
plot_bargarph_3group_comparison(pct_lossnotchosen,'pct-lossnotchosen','session1','percentage chosen',[0.4,0.7],blkname,figdir)
plot_bargarph_3group_comparison(pct_winvsloss,'pct-win vs loss chosen','session1','percentage chosen','auto',blkname,figdir)
%%%stats
 %oneway anova
 %win
 blk=2;
 [p,tbl,stats]=kruskalwallis([pct_winchosen.session1.G1.all(:,blk);...
     pct_winchosen.session1.G2.all(:,blk);...
     pct_winchosen.session1.G3.all(:,blk)],...
     groups);
 
  blk=2;
 [p,tbl,stats]=anova1([pct_winchosen.session1.G1.all(:,blk);...
     pct_winchosen.session1.G2.all(:,blk);...
     pct_winchosen.session1.G3.all(:,blk)],...
     groups);
 %post-hoc ttest
 blk=2;
 [h,p,ci,stats]=ttest2(pct_winchosen.session1.G3.all(:,blk),pct_winchosen.session1.G1.all(:,blk))
 
 %loss
  blk=3;
 [p,tbl,stats]=kruskalwallis([pct_lossnotchosen.session1.G1.all(:,blk);...
     pct_lossnotchosen.session1.G2.all(:,blk);...
     pct_lossnotchosen.session1.G3.all(:,blk)],...
     group);
 
  blk=3;
 [p,tbl,stats]=anova1([pct_lossnotchosen.session1.G1.all(:,blk);...
     pct_lossnotchosen.session1.G2.all(:,blk);...
     pct_lossnotchosen.session1.G3.all(:,blk)],...
     groups);
  %post-hoc ttest
 blk=3;
 [h,p,ci,stats]=ttest2(pct_lossnotchosen.session1.G1.all(:,blk),pct_lossnotchosen.session1.G2.all(:,blk))
 %
within=table({'1';'2';'3'},...
     'VariableNames',{'blktype'});
rm=fitrm(T,'winvslosschosen_bothv-winvslosschosen_lossv ~ group','WithinDesign',within);
ranovatbl=ranova(rm,'WithinModel','blktype');
%% plot bargarph for switch probabilities
for session={'session1'}%,'session2','session3'}
    for subset={'G1','G2','G3','all'}
        for trainset={'pos','con','all'}
            for i=1:size(switchprob.(session{1}).(subset{1}).(trainset{1}),1)
                for j=1:size(switchprob.(session{1}).(subset{1}).(trainset{1}),2)
                      switchprob_both.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_Both');
                      switchprob_nothing.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_Nothing');
                      switchprob_win.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_win');
                      switchprob_nowin.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_nowin');
                      switchprob_loss.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_loss');
                      switchprob_noloss.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_noloss');
                      switchprob_winalone.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_Winalone');
                      switchprob_lossalone.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_Lossalone');
                      switchprob_same.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_same');
                      switchprob_all.(session{1}).(subset{1}).(trainset{1})(i,j)=getfield(switchprob,session{1},subset{1},trainset{1},{i,j},'switchprob_all');
                end
            end
        end
    end
end
plot_bargarph_3group_comparison(switchprob_win,'switchprob after win','session1','switch probaility','auto',blkname,figdir)
plot_bargarph_3group_comparison(switchprob_nowin,'switchprob after nowin','session1','switch probaility','auto',blkname,figdir)
plot_bargarph_3group_comparison(switchprob_loss,'switchprob after loss','session1','switch probaility','auto',blkname,figdir)
plot_bargarph_3group_comparison(switchprob_noloss,'switchprob after noloss','session1','switch probaility','auto',blkname,figdir)
%%
plot_linegarph_3group_comparison_training_questionaires(QIDS,'QIDS',figdir)
plot_linegarph_3group_comparison_training_questionaires(tSTAI,'tSTAI',figdir)
plot_linegarph_3group_comparison_training_questionaires(sSTAI,'sSTAI',figdir)
plot_linegarph_3group_comparison_training_questionaires(RRS,'RRS',figdir)
plot_linegarph_3group_comparison_training_questionaires(DAQ,'DAQ',figdir)
% %%
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_nowin,'switchprob-nowin',figdir)
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_win,'switchprob-win',figdir)
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_loss,'switchprob-loss',figdir)
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_noloss,'switchprob-noloss',figdir)
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_nothing,'switchprob-nothing',figdir)
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_both,'switchprob-both',figdir)

%% calculate winPN loss PN and winvslossPN
nfield1=fieldnames(switchprob_win);
for field1=nfield1'
    nfield2=fieldnames(switchprob_win.(field1{1}));
    for field2=nfield2'
        nfield3=fieldnames(switchprob_win.(field1{1}).(field2{1}));
        for field3=nfield3'
        switchprob_winPN.(field1{1}).(field2{1}).(field3{1})=1-switchprob_win.(field1{1}).(field2{1}).(field3{1})-switchprob_nowin.(field1{1}).(field2{1}).(field3{1});
        switchprob_lossPN.(field1{1}).(field2{1}).(field3{1})=1-switchprob_noloss.(field1{1}).(field2{1}).(field3{1})-switchprob_loss.(field1{1}).(field2{1}).(field3{1});
        switchprob_winvslossPN.(field1{1}).(field2{1}).(field3{1})=switchprob_winPN.(field1{1}).(field2{1}).(field3{1})-switchprob_lossPN.(field1{1}).(field2{1}).(field3{1});
        end
    end
end
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_winPN,'switchprob-winPN',figdir)
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_lossPN,'switchprob-lossPN',figdir)
% plot_linegarph_3group_comparison_training_switchprobs(switchprob_winvslossPN,'switchprob-winvslossPN',figdir)

% %%
% for session={'session1'}%,'session2','session3'}
%      subset={'all'};
%       plot_bargarph_percentages_chosen_win_noloss
%       plot_bargarph_percentages_of_choosing_left_shape1
%         
%      for subset={'G1','G2','G3','all'}
%         plot_bargarph_percentages_PN
%      end
% end


% %plot bargarph with difference bet. sessions
% sessiona={'session3'};sessionb={'session2'};training={'con'};subset={'G3'};
% plot_bargarph_percentages_PN_diff_sessions
% sessiona={'session2'};sessionb={'session1'};training={'con'};subset={'G3'};
% plot_bargarph_percentages_neg_diff_sessions
% sessiona={'session2'};sessionb={'session1'};training={'con'};subset={'G3'};
% plot_bargarph_percentages_pos_diff_sessions
% corrplot(switchprob_nowin,2,tSTAI,'G1','session1');


%%
for session={'session1'}%,'session2','session3'}
allSwitch.(session{1})=asin(sqrt(switchprob_all.(session{1}).all.all));
T_all.(session{1})=array2table(allSwitch.(session{1}),'VariableNames',{['switchprob_all_bothv_',session{1}],['switchprob_all_winv_',session{1}],['switchprob_all_lossv_',session{1}]});
    
    
nowinSwitch.(session{1})=asin(sqrt(switchprob_nowin.(session{1}).all.all));
T_nowin.(session{1})=array2table(nowinSwitch.(session{1}),'VariableNames',{['switchprob_nowin_bothv_',session{1}],['switchprob_nowin_winv_',session{1}],['switchprob_nowin_lossv_',session{1}]});

winSwitch.(session{1})=asin(sqrt(switchprob_win.(session{1}).all.all));
T_win.(session{1})=array2table(winSwitch.(session{1}),'VariableNames',{['switchprob_win_bothv_',session{1}],['switchprob_win_winv_',session{1}],['switchprob_win_lossv_',session{1}]});

lossSwitch.(session{1})=asin(sqrt(switchprob_loss.(session{1}).all.all));
T_loss.(session{1})=array2table(lossSwitch.(session{1}),'VariableNames',{['switchprob0_loss_bothv_',session{1}],['switchprob_loss_winv_',session{1}],['switchprob_loss_lossv_',session{1}]});

nolossSwitch.(session{1})=asin(sqrt(switchprob_noloss.(session{1}).all.all));
T_noloss.(session{1})=array2table(nolossSwitch.(session{1}),'VariableNames',{['switchprob_noloss_bothv_',session{1}],['switchprob_noloss_winv_',session{1}],['switchprob_noloss_lossv_',session{1}]});

winPNSwitch.(session{1})=asin(sqrt(1-switchprob_win.(session{1}).all.all))-asin(sqrt(switchprob_nowin.(session{1}).all.all));
T_winPN.(session{1})=array2table(winPNSwitch.(session{1}),'VariableNames',{['switchprob_winPN_bothv_',session{1}],['switchprob_winPN_winv_',session{1}],['switchprob_winPN_lossv_',session{1}]});

lossPNSwitch.(session{1})=asin(sqrt(1-switchprob_noloss.(session{1}).all.all))-asin(sqrt(switchprob_loss.(session{1}).all.all));
T_lossPN.(session{1})=array2table(lossPNSwitch.(session{1}),'VariableNames',{['switchprob_lossPN_bothv_',session{1}],['switchprob_lossPN_winv_',session{1}],['switchprob_lossPN_lossv_',session{1}]});

winvslossPNSwitch.(session{1})=winPNSwitch.(session{1})-lossPNSwitch.(session{1});
T_winvslossPN.(session{1})=array2table(winvslossPNSwitch.(session{1}),'VariableNames',{['switchprob_winvslossPN_bothv_',session{1}],['switchprob_winvslossPN_winv_',session{1}],['switchprob_winvslossPN_lossv_',session{1}]});
histogram_3groups(QIDS,'QIDS','QIDS_histogram',figdir)
T_ques.(session{1})=array2table([tSTAI.(session{1}).all.all,sSTAI.(session{1}).all.all,QIDS.(session{1}).all.all,RRS.(session{1}).all.all,DAQ.(session{1}).all.all],...
                          'VariableNames',{['tSTAI_',session{1}],['sSTAI_',session{1}],['QIDS_',session{1}],['RRS_',session{1}],['DAQ_',session{1}]});

end
T_s=array2table(sublist.all','VariableNames',{'subject'});
T_group=array2table(group,'VariableNames',{'group'});
T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});

% T=[T_s,T_group,T_treatment,...
%     T_all.session1,T_all.session2,T_all.session3,...
%     T_nowin.session1,T_nowin.session2,T_nowin.session3,...
%     T_loss.session1,T_loss.session2,T_loss.session3,...
%     T_win.session1,T_win.session2,T_win.session3,...
%     T_noloss.session1,T_noloss.session2,T_noloss.session3,...
%     T_winPN.session1,T_winPN.session2,T_winPN.session3,...
%     T_lossPN.session1,T_lossPN.session2,T_lossPN.session3,...
%     T_winvslossPN.session1,T_winvslossPN.session2,T_winvslossPN.session3,...
%     T_ques.session1,T_ques.session2,T_ques.session3];
% writetable(T,[datadir,'swtichprobs_asinsqrt'])
% M=table(switchprob_lossPN.session1.G1.all(:,2),QIDS.session1.G1.all-QIDS.session3.G1.all);
% M.Properties.VariableNames={'switchprob_lossPN','sSTAI_inprovement'};
% writetable(M,'spss.txt','Delimiter',' ')
%%
scatterplot_3groups(sSTAI.session1.all.all,mean(switchprob_lossPN.session1.all.all,2),'state anxiety','punishment positive bias','lossPN_sSTAI','Spearman',group,figdir);
scatterplot_3groups(sSTAI.session1.all.all,mean(switchprob_winPN.session1.all.all,2),'state anxiety','reward positive bias','winPN_sSTAI','Spearman',group,figdir);
scatterplot_3groups(QIDS.session1.all.all,mean(switchprob_lossPN.session1.all.all,2),'overall depression score','punishment positive bias','lossPN_QIDS','Spearman',group,figdir);
scatterplot_3groups(QIDS.session1.all.all,mean(switchprob_winPN.session1.all.all,2),'overall depression score','reward positive bias','winPN_QIDS','Spearman',group,figdir);
% scatterplot_3groups(sSTAI.session1.all.all-sSTAI.session3.all.all,mean(switchprob_winPN.session1.all.all-switchprob_winPN.session3.all.all,2),'state anxiety session1-session3','winPN session1-session3','winPN_sSTAI_session1vssession3','Spearman',group,figdir);
% scatterplot_3groups(sSTAI.session1.all.all-sSTAI.session3.all.all,mean(switchprob_lossPN.session1.all.all-switchprob_lossPN.session3.all.all,2),'state anxiety session1-session3','lossPN session1-session3','lossPN_sSTAI_session1vssession3','Spearman',group,figdir);
%%
scatterplot_3groups(QIDS.session1.all.all,mean(switchprob_lossPN.session1.all.all,2),'overall depression score','punishment positive bias','lossPN_QIDS','Spearman',group,figdir);
scatterplot_3groups(QIDS.session1.all.all,mean(switchprob_winPN.session1.all.all,2),'overall depression score','reward positive bias','winPN_QIDS','Spearman',group,figdir);


% scatterplot_3groups(QIDS.session1.all.all-QIDS.session3.all.all,mean(switchprob_winPN.session1.all.all-switchprob_winPN.session3.all.all,2),'QIDS session1-session3','winPN session1-session3','winPN_QIDS_session1vssession3','Spearman',figdir);
% scatterplot_3groups(QIDS.session1.all.all-QIDS.session3.all.all,mean(switchprob_lossPN.session1.all.all-switchprob_lossPN.session3.all.all,2),'QIDS session1-session3','lossPN session1-session3','lossPN_QIDS_session1vssession3','Spearman',figdir);
% scatterplot_3groups(QIDS.session1.all.all,pct_winchosen.session1.all.all(:,2),'QIDS','winchosen','winchosen_QIDS','Spearman',figdir);
% scatterplot_3groups(QIDS.session1.all.all,pct_lossnotchosen.session1.all.all(:,3),'QIDS','loss not chosen','lossnotchosen_QIDS','Spearman',figdir);
% scatterplot_3groups(QIDS.session1.all.all,pct_winchosen.session1.all.all(:,2),'sSTAI','winchosen','winchosen_sSTAI','Spearman',figdir);
% scatterplot_3groups(QIDS.session1.all.all,pct_lossnotchosen.session1.all.all(:,3),'sSTAI','loss not chosen','lossnotchosen_sSTAI','Spearman',figdir);
% %%
% histogram_3groups(QIDS,'overall depression score','QIDS_histogram',figdir)
% histogram_3groups(sSTAI,'State Anxiety','sSTAI_hisrogram',figdir)
% %%
%  plot_bargarph_3group_comparison_pn(switchprob_winPN,'reward positive bias','session1',[0,0.8],blkname,figdir)
%  plot_bargarph_3group_comparison_pn(switchprob_lossPN,'punishment positive bias','session1',[0,0.8],blkname,figdir)
%  plot_bargarph_3group_comparison(switchprob_nowin,'nowin','session1',[0,0.6],blkname,figdir)
%  plot_bargarph_3group_comparison(switchprob_loss,'loss','session1',[0,0.5],blkname,figdir)
%  plot_bargarph_3group_comparison_pos(switchprob_win,'win','session1',[0.6,1],blkname,figdir)
%  plot_bargarph_3group_comparison_pos(switchprob_noloss,'noloss','session1',[0.6,1],blkname,figdir)
% 
 %%
for session={'session1'}%,'session2','session3'}
     for subset={'G1','G2','G3','all'}
     mean_switchprob_win_PSNS_diff.(session{1}).(subset{1}).all=asin(sqrt(switchprob_nowin.(session{1}).(subset{1}).all(:,2)))...
         -asin(sqrt(switchprob_win.(session{1}).(subset{1}).all(:,2)))-asin(sqrt(switchprob_nowin.(session{1}).(subset{1}).all(:,3)))...
         +asin(sqrt(switchprob_win.(session{1}).(subset{1}).all(:,3)));
     switchprob_win_PSNS_winv.(session{1}).(subset{1}).all=asin(sqrt(switchprob_nowin.(session{1}).(subset{1}).all(:,2)))...
         -asin(sqrt(switchprob_win.(session{1}).(subset{1}).all(:,2)));
     switchprob_win_PSNS_lossv.(session{1}).(subset{1}).all=asin(sqrt(switchprob_nowin.(session{1}).(subset{1}).all(:,3)))...
         -asin(sqrt(switchprob_win.(session{1}).(subset{1}).all(:,3)));
     switchprob_loss_PSNS_winv.(session{1}).(subset{1}).all=asin(sqrt(switchprob_loss.(session{1}).(subset{1}).all(:,2)))...
         -asin(sqrt(switchprob_noloss.(session{1}).(subset{1}).all(:,2)));
     switchprob_loss_PSNS_lossv.(session{1}).(subset{1}).all=asin(sqrt(switchprob_loss.(session{1}).(subset{1}).all(:,3)))...
         -asin(sqrt(switchprob_noloss.(session{1}).(subset{1}).all(:,3)));
     mean_switchprob_loss_PSNS_diff.(session{1}).(subset{1}).all=asin(sqrt(switchprob_loss.(session{1}).(subset{1}).all(:,3)))...
         -asin(sqrt(switchprob_noloss.(session{1}).(subset{1}).all(:,3)))-asin(sqrt(switchprob_loss.(session{1}).(subset{1}).all(:,2)))...
         +asin(sqrt(switchprob_noloss.(session{1}).(subset{1}).all(:,2)));
     mean_switchprob_win.(session{1}).(subset{1}).all=mean(asin(sqrt(switchprob_win.(session{1}).(subset{1}).all(:,2:3))),2);
     mean_switchprob_nowin.(session{1}).(subset{1}).all=mean(asin(sqrt(switchprob_nowin.(session{1}).(subset{1}).all(:,2:3))),2);
     mean_switchprob_loss.(session{1}).(subset{1}).all=mean(asin(sqrt(switchprob_loss.(session{1}).(subset{1}).all(:,2:3))),2);
     mean_switchprob_noloss.(session{1}).(subset{1}).all=mean(asin(sqrt(switchprob_noloss.(session{1}).(subset{1}).all(:,2:3))),2);
     
     mean_pct_win_diff.(session{1}).(subset{1}).all=asin(sqrt(pct_winchosen.(session{1}).(subset{1}).all(:,2)))-asin(sqrt(pct_winchosen.(session{1}).(subset{1}).all(:,3)));
      mean_pct_loss_diff.(session{1}).(subset{1}).all=asin(sqrt(pct_lossnotchosen.(session{1}).(subset{1}).all(:,3)))-asin(sqrt(pct_lossnotchosen.(session{1}).(subset{1}).all(:,2)));
      
     % mean_win_psns_winv_lossv=
     end
end

 plot_bargarph_3group_comparison_mean(mean_switchprob_win_PSNS_diff,'session1','Win Ppsns volatile vs stable','postive stay negative switch','auto',figdir)
  plot_bargarph_3group_comparison_mean(mean_switchprob_loss_PSNS_diff,'session1','Loss Ppsns volatile vs stable','postive stay negative switch','auto',figdir)
 plot_bargarph_3group_comparison_mean(mean_switchprob_lossPN,'session1','punishment positive bias','postive stay vs negative switch',[0,0.8],figdir)
 plot_bargarph_3group_comparison_mean(mean_switchprob_nowin,'session1','nowin','Switch probability','auto',figdir)
 plot_bargarph_3group_comparison_mean(mean_switchprob_loss,'session1','loss','Switch probability','auto',figdir)
 plot_bargarph_3group_comparison_pos_mean(mean_switchprob_win,'win','session1','auto',figdir)
 plot_bargarph_3group_comparison_mean(mean_switchprob_noloss,'session1','noloss','Switch probability','auto',figdir)
  plot_bargarph_3group_comparison_mean(mean_switchprob_win,'session1','win','Switch probability','auto',figdir)

 plot_bargarph_3group_comparison_mean(mean_pct_win_diff,'session1','Win-Chosen volatile vs stable','percentage of win chosen','auto',figdir)
  plot_bargarph_3group_comparison_mean(mean_pct_loss_diff,'session1','Loss-notchsen volatile vs stable','percentage of loss not chosen','auto',figdir)

  
   plot_bargarph_3group_comparison_mean(switchprob_win_PSNS_winv,'session1','Win Ppsns - Win Volatile','postive stay negative switch',[0 0.6],figdir)
   plot_bargarph_3group_comparison_mean(switchprob_win_PSNS_lossv,'session1','Win Ppsns - Win Stable','postive stay negative switch',[0 0.6],figdir)
   plot_bargarph_3group_comparison_mean(switchprob_loss_PSNS_winv,'session1','Loss Ppsns - Loss Stable','postive stay negative switch',[0 0.6],figdir)
   plot_bargarph_3group_comparison_mean(switchprob_loss_PSNS_lossv,'session1','Loss Ppsns - Loss Volatile','postive stay negative switch',[0 0.6],figdir)
 
   plot_bargarph_3group_comparison_mean(switchprob_loss_PSNS_lossv,'session1','Loss Ppsns - Loss Volatile','postive stay negative switch',[0 0.6],figdir)

  %  %%%stats
%  %oneway anova
%  [p,tbl,stats]=kruskalwallis([mean_pct_winvsloss.session1.G1.all,mean_pct_winvsloss.session1.G2.all,mean_pct_winvsloss.session1.G3.all]);
%  %post-hoc wilconxon rank sum test
%  [p,h,stats]=ranksum(mean_pct_winvsloss.session1.G2.all,mean_pct_winvsloss.session1.G3.all)
%  %% anova1
%  [p,tbl,stats]=kruskalwallis([mean_switchprob_winPN.session1.G1.all,mean_switchprob_winPN.session1.G2.all,mean_switchprob_winPN.session1.G3.all])
%  [p,tbl,stats]=kruskalwallis([mean_switchprob_lossPN.session1.G1.all,mean_switchprob_lossPN.session1.G2.all,mean_switchprob_lossPN.session1.G3.all])
%  [p,tbl,stats]=anova1([mean_switchprob_noloss.session1.G1.all;mean_switchprob_noloss.session1.G2.all;mean_switchprob_noloss.session1.G3.all],groups)
% [p,tbl,stats]=anova1([mean_switchprob_loss.session1.G1.all;mean_switchprob_loss.session1.G2.all;mean_switchprob_loss.session1.G3.all],groups)
% [p,tbl,stats]=anova1([mean_switchprob_nowin.session1.G1.all;mean_switchprob_nowin.session1.G2.all;mean_switchprob_nowin.session1.G3.all],groups)
% [p,tbl,stats]=anova1([mean_switchprob_win.session1.G1.all;mean_switchprob_win.session1.G2.all;mean_switchprob_win.session1.G3.all],groups)
% 
% scatterplot_3groups(tSTAI.session1.all.all,mean_switchprob_win.session1.all.all,'-','-','noloss_switch','Spearman',groups,figdir);
% scatterplot_3groups(sSTAI.session1.all.all,mean_switchprob_win.session1.all.all,'-','-','noloss_switch','Spearman',groups,figdir);
% scatterplot_3groups(DAQ.session1.all.all,mean_switchprob_win.session1.all.all,'-','-','noloss_switch','Spearman',groups,figdir);
% 
% scatterplot_3groups(tSTAI.session1.all.all,mean_switchprob_noloss.session1.all.all,'-','-','noloss_switch','Pearson',groups,figdir);
% scatterplot_3groups(sSTAI.session1.all.all,mean_switchprob_noloss.session1.all.all,'-','-','noloss_switch','Pearson',groups,figdir);
% 
% %%
% T_s=array2table(sublist.all','VariableNames',{'subject'});
% T_group=array2table(groups,'VariableNames',{'group'});
% T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});
% T_switch=array2table([asin(sqrt(switchprob_win.session1.all.all(:,2:3))),asin(sqrt(switchprob_nowin.session1.all.all(:,2:3))),...
%     asin(sqrt(switchprob_noloss.session1.all.all(:,2:3))),asin(sqrt(switchprob_loss.session1.all.all(:,2:3)))],...
%     'VariableNames',{'win_switch_winv','win_switch_lossv',...
%                                 'nowin_switch_winv','nowin_switch_lossv',...
%                                 'noloss_switch_winv','noloss_switch_lossv',...
%                                 'loss_switch_winv','loss_switch_lossv'});
% T=[T_s,T_group,T_treatment,T_switch];
% % within=table(repelem({'win';'loss'},3,1),repmat({'winv';'winv';'wins'},2,1),repmat({'lossv';'losss';'lossv'},2,1),...
% %      'VariableNames',{'valence','winvol','lossvol'});
% within=table(repelem({'win';'loss'},4,1),repmat({'pos';'pos';'neg';'neg'},2,1),repmat({'winv';'lossv'},4,1),...
%      'VariableNames',{'valence','outcometype','blktype'});
% rm=fitrm(T,'win_switch_winv-loss_switch_lossv ~ group','WithinDesign',within);
% ranovatbl=ranova(rm,'WithinModel','valence*outcometype*blktype');
% %%
% T_s=array2table(sublist.all','VariableNames',{'subject'});
% T_group=array2table(groups,'VariableNames',{'group'});
% T_treatment=array2table(trainingtype.all,'VariableNames',{'treatment'});
% T_switch=array2table([mean_switchprob_win.session1.all.all,mean_switchprob_nowin.session1.all.all,...
%     mean_switchprob_noloss.session1.all.all,mean_switchprob_loss.session1.all.all],...
%     'VariableNames',{'win_switch',...
%                                 'nowin_switch',...
%                                 'noloss_switch',...
%                                 'loss_switch'});
% T=[T_s,T_group,T_treatment,T_switch];
% % within=table(repelem({'win';'loss'},3,1),repmat({'winv';'winv';'wins'},2,1),repmat({'lossv';'losss';'lossv'},2,1),...
% %      'VariableNames',{'valence','winvol','lossvol'});
% within=table(repelem({'win';'loss'},2,1),repmat({'pos';'neg'},2,1),...
%      'VariableNames',{'valence','outcometype'});
% rm=fitrm(T,'win_switch-loss_switch ~ group','WithinDesign',within);
% ranovatbl=ranova(rm,'WithinModel','valence*outcometype');
% 
scatterplot_3groups(loss_adpt_vol_adpt_model.session1.all.all,switchprob_noloss.session1.all.all(:,3)-switchprob_noloss.session1.all.all(:,1),'loss learning rate adptation','no loss vol vs stable',...
    'scatterplot_loss_adpt_corr_switch_noloss','Pearson',groups,figdir)
scatterplot_3groups(loss_adpt_vol_adpt_model.session1.all.all,switchprob_loss.session1.all.all(:,3)-switchprob_loss.session1.all.all(:,1),'loss learning rate adptation','loss vol vs stable',...
    'scatterplot_loss_adpt_corr_switch_loss','Pearson',groups,figdir)
scatterplot_3groups(win_adpt_vol_adpt_model.session1.all.all,switchprob_win.session1.all.all(:,3)-switchprob_win.session1.all.all(:,1),'win learning rate adptation','win vol vs stable',...
    'scatterplot_win_adpt_corr_switch_win','Pearson',groups,figdir)
scatterplot_3groups(win_adpt_vol_adpt_model.session1.all.all,switchprob_nowin.session1.all.all(:,3)-switchprob_nowin.session1.all.all(:,1),'win learning rate adptation','nowin vol vs stable',...
    'scatterplot_win_adpt_corr_switch_nowin','Pearson',groups,figdir)

