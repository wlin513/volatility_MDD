function plot_pupil_for_each_block(inpwin,inploss,session,inptype,blkname,figdir)
    %%this function plot win and loss pupil responses seperately for each block for the clinical study
    %the inputs are windata,lossdata,session name,input type(e.g. received vs not, received, not
    %received),blkname,save destination
    
    % plot winvol winstable lossvol lossstable across blocks\
    if nargin<5
        blkname={'both volatile','win volatile','loss volatile'};
    end
    if nargin<6
        figdir= '/home/wlin/Documents/2017_clinical_study/tmp_fig/';
    end
    atr.color.G1=[0.5 0 0.8];
    atr.color.G2=[0.19 0.73 0.74];
    atr.color.G3=[1.0 0.6 0.2];
    atr.groupname.G1='depression group';
    atr.groupname.G2='remitted depression group';
    atr.groupname.G3='control group';
    f1=figure('units','inch','position',[0,0,16,8]);
    for blk=1:3

        xpoints=1:3500;
        y(1,1:length(xpoints))=zeros; %%reference line
        %% plot pupil response for win delta
        subplot(2,3,blk)
        for group={'G1','G2','G3'}
            means=mean(inpwin.(session).(group{1}).all(:,blk,:));
            stes=std(inpwin.(session).(group{1}).all(:,blk,:))./sqrt(size(inpwin.(session).(group{1}).all(:,blk,:),1));

            upper=reshape(means(:,:,501:end)+stes(:,:,501:end),1,3500);
            lower=reshape(means(:,:,501:end)-stes(:,:,501:end),1,3500);
            jbfill(xpoints,upper,lower,atr.color.(group{1}),atr.color.(group{1}),0, 0.15);
            hold on
            plot(xpoints,reshape(means(:,:,501:end),1,3500),'Color',atr.color.(group{1}),'LineWidth',2,'DisplayName',atr.groupname.(group{1}));
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
        end
           %
            ax=gca;
            ax.XTick=[0:500:3500];
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000','6000'};
            xlim([0 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel(['Pupil response - ',inptype]);
            legend('Location',[-0.06+0.285*blk,0.85,0.14,0.06])
            title(['win outcomes in ',blkname{blk},' block'])
            hold on
        %%
            subplot(2,3,blk+3);
        for group={'G1','G2','G3'}
            means=mean(inploss.(session).(group{1}).all(:,blk,:));
            stes=std(inploss.(session).(group{1}).all(:,blk,:))./sqrt(size(inpwin.(session).(group{1}).all(:,blk,:),1));

            upper=reshape(means(:,:,501:end)+stes(:,:,501:end),1,3500);
            lower=reshape(means(:,:,501:end)-stes(:,:,501:end),1,3500);
            jbfill(xpoints,upper,lower,atr.color.(group{1}),atr.color.(group{1}),0, 0.15);
            hold on
            plot(xpoints,reshape(means(:,:,501:end),1,3500),'Color',atr.color.(group{1}),'LineWidth',2,'DisplayName',atr.groupname.(group{1}));
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
        end
               %
            ax=gca;
            ax.XTick=[0:500:3500];
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000','6000'};
            xlim([0 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel(['Pupil response - ',inptype]);
            legend('Location',[-0.06+0.285*blk,0.37,0.14,0.06])
            title(['loss outcomes in ',blkname{blk},' block'])
    end
   saveas(f1,[figdir,'pupil/',session,'_win_loss_',inptype,'_for_each_block.png']);