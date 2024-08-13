%% 
clear all;
clc;
%%
cd('/home/wlin/Documents/2017_clinical_study/code/pupil/elife_analysis');
getfolders;
%% load subject data pool

session1_sublist.G1=1001:1041;
session1_sublist.G2=2001:2042; session1_sublist.G2(session1_sublist.G2==2033)=[];
session1_sublist.G3=3001:3041; 
% exclusion: exclude participant who some of the data were missing
%G1
session1_sublist.G1(session1_sublist.G1==1005)=[];%no session 1 pupil data
%session1_sublist.G1(session1_sublist.G1==1009)=[];%no session 2 pupil data
%session1_sublist.G1(session1_sublist.G1==1014)=[];%no session 3 pupil data

%G2
%session1_sublist.G2(session1_sublist.G2==2012)=[];%no session 2 pupil data
%session1_sublist.G2(session1_sublist.G2==2016)=[];%no session 2 pupil data
%G3
session1_sublist.G3(session1_sublist.G3==3040)=[];%couldn't read event asc file but edf2mat file is fine. double check when have time/needed
session1_sublist.all=[session1_sublist.G1,session1_sublist.G2,session1_sublist.G3];
%%

for sub=1:length(session1_sublist.all)
    sub
    tracker_chop_session1_byblock(num2str(session1_sublist.all(sub)),datadir);
end