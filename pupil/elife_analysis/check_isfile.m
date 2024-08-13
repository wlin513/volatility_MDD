

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
sublist.all=[sublist.G1,sublist.G2,sublist.G3];
%%
notcompletedatalist.s1=nan;
notcompletedatalist.s2=nan;
notcompletedatalist.s3=nan;
ns1=1;ns2=1;ns3=1;
for sub=1:size(sublist.all,2)
    if ~isfile([datadir,num2str(sublist.all(sub)),'/',num2str(sublist.all(sub)),'_session_1_vol_linked_sample.asc'])
        notcompletedatalist.s1(ns1)=sublist.all(sub);
        ns1=ns1+1;
    end
    if ~isfile([datadir,num2str(sublist.all(sub)),'/',num2str(sublist.all(sub)),'_session_2_vol_linked_sample.asc'])
        notcompletedatalist.s2(ns2)=sublist.all(sub);
        ns2=ns2+1;
    end
    if ~isfile([datadir,num2str(sublist.all(sub)),'/',num2str(sublist.all(sub)),'_session_3_vol_linked_sample.asc'])
        notcompletedatalist.s3(ns3)=sublist.all(sub);
        ns3=ns3+1;
    end
end