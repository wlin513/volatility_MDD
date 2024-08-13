%%check the eye tracking data for missing parts
clear all
getfolders
%% load subject data pool
sublist.G1=1001:1041;sublist.G1(sublist.G1==1036)=[];
sublist.G2=2001:2042; sublist.G2(sublist.G2==2033)=[];sublist.G2(sublist.G2==2019)=[];
sublist.G3=3001:3041; sublist.G3(sublist.G3==3033)=[];
sublist.G1(sublist.G1==1004)=[];%exclueded from beh analysis
sublist.G1(sublist.G1==1039)=[];%exclueded from beh analysis
sublist.G2(sublist.G1==2034)=[];%exclueded from beh analysis
% exclude participants for session 1
sublist.G1(sublist.G1==1005)=[];%no session 1 pupil data
%interpolated more than 60% of the pupil data
sublist.G1(sublist.G1==1031)=[];% s1 62% s3 72%
sublist.G3(sublist.G3==3024)=[];%somehow all session1 results are nans
sublist.G3(sublist.G3==3009)=[];%s1 65% s3 65% 
sublist.G3(sublist.G3==3029)=[];%s1 72% s3 75% 
sublist.G3(sublist.G3==3040)=[];%couldn't read event asc file but edf2mat file is fine. double check when have time/needed
sublist.all=[sublist.G1,sublist.G2,sublist.G3];
group=[ones(size(sublist.G1,2),1);ones(size(sublist.G2,2),1)*2;ones(size(sublist.G3,2),1)*3];
%%
sublist.G1=1036;
sublist.G2=2019;
sublist.G3=3033;
sublist.all=[sublist.G1,sublist.G2,sublist.G3];
%%
for sub=1:size(sublist.all,2)
        for session={'session_1'}%,'session_2','session_3'}
            name=num2str(sublist.all(sub));
            file_name=['out_',name,'_',session{1},'_normalized_within_session1.mat'];
            load([datadir,name,'/',file_name])

            sub_data.tracker_data=out; 
            store=parsemisstracker(sub_data);


            for i=1:size(store.options_missa,1)    
                rewards_missing(sub,i)=sum(store.rew_outcome_missb(i,:))/size(store.options_missa,2)*100;
                loss_missing(sub,i)=sum(store.pun_outcome_missb(i,:))/size(store.options_missa,2)*100;   
            end

            rew_exclude(1,240)=zeros;
            loss_exclude(1,240)=zeros;
            for k=1:size(store.options_missa,1) 

                if rewards_missing(sub,k)<=50
                    rew_exclude(k)=0;
                else
                    rew_exclude(k)=1;
                end

                if loss_missing(sub,k)<=50
                    pun_exclude(k)=0;
                else
                    pun_exclude(k)=1;
                end
            end

            save([datadir,name,'/',name,'_',session{1}, '_pupil_exclude_trials.mat'],'rew_exclude','pun_exclude')
        end
end

% clear all
% subs={'v9','v10','v11','v12','v13','v15','v16','v18', 'v21','v23','v25','v26', 'v28','v29','v31','v34','v38','v40','v42','v43','v44','v45','v47','v48','v51','v52','v54','v55','v56'};
% for sub=1:size(subs,2)
%     
% name=subs{sub};
% dir_name=['C:\Users\epulcu\Desktop\Oxford_Volatility_MDD\',name,'\'];
% 
% pth=dir_name;
% cd(pth)
% load out.mat
% 
% pth='C:\Users\epulcu\Desktop\Oxford_Volatility_MDD\';
% cd(pth)
% 
% sub_data.tracker_data=out;
% store=parsemisstracker(sub_data);
% 
% 
% for i=1:size(store.options_missa,1)    
%     rewards_missing(sub,i)=sum(store.rew_outcome_missb(i,:))/size(store.options_missa,2)*100;
%     loss_missing(sub,i)=sum(store.pun_outcome_missb(i,:))/size(store.options_missa,2)*100;   
% end
% 
% 
% rew_exclude=[];
% loss_exclude=[];
% for k=1:size(store.options_missa,1) 
%     
%     if rewards_missing(sub,k)<=50
%         rew_exclude(k)=0;
%     else
%         rew_exclude(k)=1;
%     end
%     
%     if loss_missing(sub,k)<=50
%         pun_exclude(k)=0;
%     else
%         pun_exclude(k)=1;
%     end
% end
% name=subs{sub};
% dir_name=['C:\Users\epulcu\Desktop\Oxford_Volatility_MDD\',name,'\'];
% 
% pth=dir_name;
% cd(pth)
% 
% save([name '_pupil_exclude_trials.mat'],'rew_exclude','pun_exclude')
% 
% end
% 
