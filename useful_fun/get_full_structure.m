function out=get_full_structure(inp,trainingtype,sublist)

 out=inp;
for session=fieldnames(inp)'
    out.(session{1}).all.pos=inp.(session{1}).all.all(trainingtype.all==1,:,:);
    out.(session{1}).all.con=inp.(session{1}).all.all(trainingtype.all==2,:,:);
     for subset={'G1','G2','G3'}
        subspool=sublist.(subset{1});
        out.(session{1}).(subset{1}).all=out.(session{1}).all.all(ismember(sublist.all,subspool),:,:);
        out.(session{1}).(subset{1}).pos=out.(session{1}).(subset{1}).all(trainingtype.(subset{1})==1,:,:);
        out.(session{1}).(subset{1}).con=out.(session{1}).(subset{1}).all(trainingtype.(subset{1})==2,:,:);
    end
end