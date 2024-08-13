function out=get_p1vital_info(dataout,sublist,trainingtype,target)
sn=1;
snames={'baseline','visit2','followup'};
for session={'session1','session2','session3'}
     subspool=sublist.all;
        for ss=1:size(subspool,2)
        subi=find(strcmp(dataout.subject_names,['GB-25-',num2str(subspool(ss))]));
             if isfield(dataout.subject_data(subi).(target),snames{sn})
                out.(session{1}).all.all(ss,1)=dataout.subject_data(subi).(target).(snames{sn});
             else 
                 out.(session{1}).all.all(ss,1)=nan;
            end
        end 
  sn=sn+1;
end
out=get_full_structure(out,trainingtype,sublist);
