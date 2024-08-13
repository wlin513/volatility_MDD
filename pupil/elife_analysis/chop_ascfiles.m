%% 
clear all;
clc;
%%
cd('/home/wlin/Documents/2017_clinical_study/code/pupil/elife_analysis');
getfolders;
%% load subject data pool

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
%%

for sub=1:length(sublist.G3)
    tracker_chop_normalise_amongst_visits(num2str(sublist.G3(sub)),datadir);
end