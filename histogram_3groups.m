function histogram_3groups(ques,ttl,picname,figdir)
fhis=figure;
barvar=[];
newbarvar=[];
for i=min(ques.session1.all.all):max(ques.session1.all.all)
    barvar=[barvar;sum(ques.session1.G1.all==i),sum(ques.session1.G2.all==i),sum(ques.session1.G3.all==i)];
end

if size(barvar,1)>30
    for i=2:10
        if mod(size(barvar,1),i)==0
            wid=i;
        end
    end
    
    for i=1:size(barvar,1)/wid
        newbarvar=[newbarvar;sum(barvar((i-1)*wid+1:i*wid,:))];
    end
    
    b=bar(newbarvar,'stacked');
   % xticks
   xt=xticks;
   xticks(xt);
   xticklabels(xt*wid);
else
    b=bar(barvar,'stacked');
end

b(1).FaceColor=[0.5 0 0.8];
b(2).FaceColor=[0.19 0.73 0.74];
b(3).FaceColor=[1.0 0.6 0.2];
for i=1:3
    b(i).EdgeColor='none';
end
legend('depression','remitted depression','control')
legend boxoff
title(['\fontsize{16}',ttl]);
print(fhis,[figdir,picname,'.png'],'-dpng','-r300')