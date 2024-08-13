function [r,p,fcorr]=corrplot(switchp,i,ques,subset,ssession)
clc
[r,p]=corr(switchp.(ssession).(subset).all(:,i),ques.session1.(subset).all-ques.session3.(subset).all)
fcorr=figure;
scatter(switchp.(ssession).(subset).all(:,i),ques.session1.(subset).all-ques.session3.(subset).all)
lsline
