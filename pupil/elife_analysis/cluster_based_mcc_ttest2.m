function [orisigpoints, bootsigpoints, clustersigpoints] = cluster_based_mcc_ttest2(data1,data2,sigalpha,b)

%two sample t-test cluster based multiple comparison correction
%
%input should be (1)data1(sub,time), (2)data2(sub,time), (3)significant alpha level, (4)bootstrapping number;
%output would be (1)original significant points(uncorrected student t test), (2)bootstrapping significant points(uncorrected bootstrapping test), 
%(3)significant pointed after cluster based multiple comparision correction


[H,orip,CI,stats] = ttest2(data1,data2);
orit = stats.tstat;

orisigpoints=(orip<=sigalpha);
orisigpoints=double(orisigpoints);
orisigpoints(orisigpoints==0)=nan;
%% do bootstrap
Nf = size(data1,2);

th = b*(1-sigalpha);
% lo = round(b*sigalpha./2);
% hi = b-lo;
% lo = lo+1;

boott = zeros(b,Nf);
bootp = boott;
cdata1 = data1-repmat(mean(data1),[size(data1,1) 1]);
cdata2 = data2-repmat(mean(data2),[size(data2,1) 1]);
for kk=1:b
    [H,bootp(kk,:),CI,stats] = ttest2(cdata1(randi(size(data1,1),size(data1,1),1),:),cdata2(randi(size(data2,1),size(data2,1),1),:));
    boott(kk,:) = stats.tstat;
end
%% get bootstrap univariate thresholds

% IMPORTANT: as you will see, on it's own the bootstrap DOES NOT
%            correct for multiple comparisons

sortboot = sort(boott.^2); % we use t.^2 (=F) to make our life easier
unibootth = sortboot(th,:); 

bootsigpoints=(orit.^2>=unibootth);
bootsigpoints=double(bootsigpoints);
bootsigpoints(bootsigpoints==0)=nan;
%% get cluster based significant points
% get cluster threshold

clusterbootth = limo_ecluster_make(boott'.^2,bootp',sigalpha);

% cluster original data

sigcluster = limo_ecluster_test(orit.^2,orip,clusterbootth,sigalpha);

% get cluster based significant points
clustersigpoints=sigcluster.elec_mask;
clustersigpoints(clustersigpoints==0)=nan;


