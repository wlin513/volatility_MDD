%%check the eye tracking data for missing parts
%%check the eye tracking data for missing parts
clear all
getfolders
%%
sublist.G1=1001:1041;sublist.G1(sublist.G1==1036)=[];
sublist.G2=2001:2042; sublist.G2(sublist.G2==2033)=[];sublist.G2(sublist.G2==2019)=[];
sublist.G3=3001:3041; sublist.G3(sublist.G3==3033)=[];
% exclusion: exclude participant who some of the data were missing
%G1
sublist.G1(sublist.G1==1005)=[];%no session 1 pupil data
sublist.G1(sublist.G1==1009)=[];%no session 2 pupil data
sublist.G1(sublist.G1==1014)=[];%no session 3 pupil data

%G2
sublist.G2(sublist.G2==2012)=[];%no session 2 pupil data
sublist.G2(sublist.G2==2016)=[];%no session 2 pupil data
%G3
sublist.G3(sublist.G3==3040)=[];%couldn't read event asc file but edf2mat file is fine. double check when have time/needed
sublist.all=[sublist.G1,sublist.G2,sublist.G3];
group=[ones(size(sublist.G1,2),1);ones(size(sublist.G2,2),1)*2;ones(size(sublist.G3,2),1)*3];

%%
for sub=1:size(sublist.all,2)
        for session={'session_1','session_2','session_3'}
            name=num2str(sublist.all(sub));
            file_name=['out_',name,'_',session{1},'_normalised_amongst_3sessions.mat'];
            load([datadir,name,'/',file_name])

            removed=out.removed_data;

            counter=sum(removed(:,1)&removed(:,2));
            missing_percent.(session{1})(sub,:)=[sublist.all(sub),(counter/size(removed,1))*100];
        end
end
% 
% %%check the eye tracking data for missing parts- mag task
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
% removed=out.removed_data;
% 
% counter=0;
% for i=1:size(removed,1)    
%     if removed(i,1)==1 && removed(i,2)==1 
%         counter=counter+1;
%     end
% end
% 
% missing_percent(sub,1)=(counter/size(removed,1))*100;
% end
% 
