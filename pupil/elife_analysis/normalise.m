function out=normalise(in)

% normalise a vector or matrix of data


mean_in=mean(in);

std_in=std(in,0,1);

out=zeros(size(in));

for i=1:size(in,2)
out(:,i)=in(:,i)-mean_in(i);
out(:,i)=out(:,i)./std_in(i);
end