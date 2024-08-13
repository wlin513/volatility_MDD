function [orisigpoints, bootsigpoints, clustersigpoints] = cluster_based_mcc_anova1(data,group,sigalpha,b)

%two sample t-test cluster based multiple comparison correction
%
%input should be (1)data1(sub,time), (2)group information, (3)significant alpha level, (4)bootstrapping number;
%output would be (1)original significant points(uncorrected student t test), (2)bootstrapping significant points(uncorrected bootstrapping test), 
%(3)significant pointed after cluster based multiple comparision correction
for i=1:size(data,2)
 [orip(i),rbl]=anova1(data(:,i),group,'off');
 orif(i)=rbl{2,5};
end

orisigpoints=(orip<=sigalpha);
orisigpoints=double(orisigpoints);
orisigpoints(orisigpoints==0)=nan;
%% do bootstrap
Nf = size(data,2);

th = b*(1-sigalpha);
% lo = round(b*sigalpha./2);
% hi = b-lo;
% lo = lo+1;

bootf = zeros(b,Nf);
bootp = bootf;
cdata = data-repmat(mean(data),[size(data,1) 1]);
for kk=1:b
    permdata=cdata(randi(size(data,1),size(data,1),1),:);
    for i=1:size(data,2)
       [bootp(kk,i),rbl]=anova1(permdata(:,i),group,'off');
       bootf(kk,i)=rbl{2,5};
    end
end
%% get bootstrap univariate thresholds

% IMPORTANT: as you will see, on it's own the bootstrap DOES NOT
%            correct for multiple comparisons

sortboot = sort(bootf); % we use t.^2 (=F) to make our life easier
unibootth = sortboot(th,:); 

bootsigpoints=(orif>=unibootth);
bootsigpoints=double(bootsigpoints);
bootsigpoints(bootsigpoints==0)=nan;
%% get cluster based significant points
% get cluster threshold

clusterbootth = limo_ecluster_make(bootf',bootp',sigalpha);

% cluster original data

sigcluster = limo_ecluster_test(orif,orip,clusterbootth,sigalpha);

% get cluster based significant points
clustersigpoints=sigcluster.elec_mask;
clustersigpoints(clustersigpoints==0)=nan;


