%% plot winvol winstable lossvol lossstable across blocks\
for group={'reboxetine','placebo'}
      delta_win_bothvvslossv.(group{1})=squeeze(delta_win.(group{1})(:,:,1,:)-delta_win.(group{1})(:,:,3,:));
      delta_win_winvvsboths.(group{1})=squeeze(delta_win.(group{1})(:,:,2,:)-delta_win.(group{1})(:,:,4,:));
      delta_loss_bothvvswinv.(group{1})=squeeze(delta_loss.(group{1})(:,:,1,:)-delta_loss.(group{1})(:,:,2,:));
      delta_loss_lossvvsboths.(group{1})=squeeze(delta_loss.(group{1})(:,:,3,:)-delta_loss.(group{1})(:,:,4,:));
      group_means_delta_win_bothvvslossv.(group{1})=squeeze(mean(delta_win_bothvvslossv.(group{1}),2));
      group_means_delta_win_winvvsboths.(group{1})=squeeze(mean(delta_win_winvvsboths.(group{1}),2));
      group_means_delta_loss_bothvvswinv.(group{1})=squeeze(mean(delta_loss_bothvvswinv.(group{1}),2));
      group_means_delta_loss_lossvvsboths.(group{1})=squeeze(mean(delta_loss_lossvvsboths.(group{1}),2));
      group_stes_delta_win_bothvvslossv.(group{1})=squeeze(std(delta_win_bothvvslossv.(group{1}),0,2)./sqrt(length(sublist.(group{1}))));
      group_stes_delta_win_winvvsboths.(group{1})=squeeze(std(delta_win_winvvsboths.(group{1}),0,2)./sqrt(length(sublist.(group{1}))));
      group_stes_delta_loss_bothvvswinv.(group{1})=squeeze(std(delta_loss_bothvvswinv.(group{1}),0,2)./sqrt(length(sublist.(group{1}))));
      group_stes_delta_loss_lossvvsboths.(group{1})=squeeze(std(delta_loss_lossvvsboths.(group{1}),0,2)./sqrt(length(sublist.(group{1}))));
      %visit2 - visit1
      diff_delta_win_bothvvslossv.(group{1})=squeeze(delta_win_bothvvslossv.(group{1})(2,:,:)-delta_win_bothvvslossv.(group{1})(1,:,:));
      diff_delta_win_winvvsboths.(group{1})=squeeze(delta_win_winvvsboths.(group{1})(2,:,:)-delta_win_winvvsboths.(group{1})(1,:,:));
      diff_delta_loss_bothvvswinv.(group{1})=squeeze(delta_loss_bothvvswinv.(group{1})(2,:,:)-delta_loss_bothvvswinv.(group{1})(1,:,:));
      diff_delta_loss_lossvvsboths.(group{1})=squeeze(delta_loss_lossvvsboths.(group{1})(2,:,:)-delta_loss_lossvvsboths.(group{1})(1,:,:));
      group_means_diff_delta_win_bothvvslossv.(group{1})=mean(diff_delta_win_bothvvslossv.(group{1}));
      group_means_diff_delta_win_winvvsboths.(group{1})=mean(diff_delta_win_winvvsboths.(group{1}));
      group_means_diff_delta_loss_bothvvswinv.(group{1})=mean(diff_delta_loss_bothvvswinv.(group{1}));
      group_means_diff_delta_loss_lossvvsboths.(group{1})=mean(diff_delta_loss_lossvvsboths.(group{1}));
      group_stes_diff_delta_win_bothvvslossv.(group{1})=squeeze(std(diff_delta_win_bothvvslossv.(group{1}))./sqrt(length(sublist.(group{1}))));
      group_stes_diff_delta_win_winvvsboths.(group{1})=squeeze(std(diff_delta_win_winvvsboths.(group{1}))./sqrt(length(sublist.(group{1}))));
      group_stes_diff_delta_loss_bothvvswinv.(group{1})=squeeze(std(diff_delta_loss_bothvvswinv.(group{1}))./sqrt(length(sublist.(group{1}))));
      group_stes_diff_delta_loss_lossvvsboths.(group{1})=squeeze(std(diff_delta_loss_lossvvsboths.(group{1}))./sqrt(length(sublist.(group{1}))));
  
      for visit=1:2

            f1=figure;%main figure combining both time courses [reaction to chosen rewards]
            xpoints=1:3000;
            y(1,1:length(xpoints))=zeros; %%reference line
           
            %% plot win pupil response 
            %  win bothv-lossv
            upper=group_means_delta_win_bothvvslossv.(group{1})(visit,501:end)+group_stes_delta_win_bothvvslossv.(group{1})(visit,501:end);
            lower=group_means_delta_win_bothvvslossv.(group{1})(visit,501:end)-group_stes_delta_win_bothvvslossv.(group{1})(visit,501:end);
            jbfill(xpoints,upper,lower,[.15 0.35 0.26],[.15 0.35 0.26],0, 0.15);
            hold on
            plot(xpoints,group_means_delta_win_bothvvslossv.(group{1})(visit,501:end),'Color',[.15 0.35 0.26],'LineWidth',2,'DisplayName','loss volatile');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off'); 
            hold on
            %win winv-boths
            upper=group_means_delta_win_winvvsboths.(group{1})(visit,501:end)+group_stes_delta_win_winvvsboths.(group{1})(visit,501:end);
            lower=group_means_delta_win_winvvsboths.(group{1})(visit,501:end)-group_stes_delta_win_winvvsboths.(group{1})(visit,501:end);
            jbfill(xpoints,upper,lower,[0.77 0.91 0.60],[0.77 0.91 0.60],0, 0.15);
            hold on
            plot(xpoints,group_means_delta_win_winvvsboths.(group{1})(visit,501:end),'Color',[0.77 0.91 0.60],'LineWidth',2,'DisplayName','loss stable');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
           %%
            ax=gca;
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000'};
            %xlim([500 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel('Pupil response to chosen vs unchosen outcomes');
%             set( gca                       , ...
%                 'FontName'   , 'Helvetica' );
%             set([ hXLabel, hYLabel], ...
%                 'FontName'   , 'Helvetica');
% 
%             set([hXLabel, hYLabel]  , ...
%                 'FontSize'   , 13          );

            set(gca, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.02 .02] , ...
                'XMinorTick'  , 'on'      , ...
                'YMinorTick'  , 'on'      , ...
                'YGrid'       , 'off'      , ...
                'XColor'      , [0 0 0 ], ...
                'YColor'      , [0 0 0], ...
                'LineWidth'   , 2         );
            %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8/1.5 6/1.5])
            axis square
            hold on
            legend
            title([group{1},': visit',num2str(visit)])
            saveas(f1,[figdir,'pupil/',group{1},'_visit',num2str(visit),'_win_pupil_response_vol_vs_sta_to_othervol.png']);
            %% -------
            f2=figure;
            %xpoints=1:length(lossvol_lossv);
            y(1,1:length(xpoints))=zeros; %%reference line

            %% plot loss pupil response 
            %  loss bothv-winv
            upper=group_means_delta_loss_bothvvswinv.(group{1})(visit,501:end)+group_stes_delta_loss_bothvvswinv.(group{1})(visit,501:end);
            lower=group_means_delta_loss_bothvvswinv.(group{1})(visit,501:end)-group_stes_delta_loss_bothvvswinv.(group{1})(visit,501:end);
            jbfill(xpoints,upper,lower,[0.31 0.15 0.27],[0.31 0.15 0.27],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_delta_loss_bothvvswinv.(group{1})(visit,501:end),1,3000),'Color',[0.31 0.15 0.27],'LineWidth',2,'DisplayName','win volatile');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
            %loss lossv-boths
            upper=group_means_delta_loss_lossvvsboths.(group{1})(visit,501:end)+group_stes_delta_loss_lossvvsboths.(group{1})(visit,501:end);
            lower=group_means_delta_loss_lossvvsboths.(group{1})(visit,501:end)-group_stes_delta_loss_lossvvsboths.(group{1})(visit,501:end);
            jbfill(xpoints,upper,lower,[1 0.6 0.60],[1 0.6 0.60],0, 0.15);
            hold on
            plot(xpoints,group_means_delta_loss_lossvvsboths.(group{1})(visit,501:end),'Color',[1 0.6 0.60],'LineWidth',2,'DisplayName','win stable');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');        
          %%
            ax=gca;
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000'};
            %xlim([500 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel('Pupil response to chosen vs unchosen outcomes');
%             set( gca                       , ...
%                 'FontName'   , 'Helvetica' );
%             set([ hXLabel, hYLabel], ...
%                 'FontName'   , 'Helvetica');
% 
%             set([hXLabel, hYLabel]  , ...
%                 'FontSize'   , 13          );

            set(gca, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.02 .02] , ...
                'XMinorTick'  , 'on'      , ...
                'YMinorTick'  , 'on'      , ...
                'YGrid'       , 'off'      , ...
                'XColor'      , [0 0 0 ], ...
                'YColor'      , [0 0 0], ...
                'LineWidth'   , 2         );
            %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8/1.5 6/1.5])
            axis square
            hold on
            %legend([f1(1),f1(2),f1(3),f1(4)],'Win Volatile','Win Stable','Loss Volatile','Loss Stable')
            legend
            title([group{1},': visit',num2str(visit)])

            saveas(f2,[figdir,'pupil/',group{1},'_visit',num2str(visit),'_loss_pupil_response_vol_vs_sta_to_othervol.png']);
      end
    %% plot win vol visit2 - visit1
            f3=figure;
            %% plot win pupil response for each block 
            %  win bothv-lossv
            upper=group_means_diff_delta_win_bothvvslossv.(group{1})(1,501:end)+group_stes_diff_delta_win_bothvvslossv.(group{1})(1,501:end);
            lower=group_means_diff_delta_win_bothvvslossv.(group{1})(1,501:end)-group_stes_diff_delta_win_bothvvslossv.(group{1})(1,501:end);
            jbfill(xpoints,upper,lower,[.15 0.35 0.26],[.15 0.35 0.26],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_delta_win_bothvvslossv.(group{1})(1,501:end),'Color',[.15 0.35 0.26],'LineWidth',2,'DisplayName','loss volatile');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
            %win winv-boths
            upper=group_means_diff_delta_win_winvvsboths.(group{1})(1,501:end)+group_stes_diff_delta_win_winvvsboths.(group{1})(1,501:end);
            lower=group_means_diff_delta_win_winvvsboths.(group{1})(1,501:end)-group_stes_diff_delta_win_winvvsboths.(group{1})(1,501:end);
            jbfill(xpoints,upper,lower,[0.77 0.91 0.60],[0.77 0.91 0.60],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_delta_win_winvvsboths.(group{1})(1,501:end),'Color',[0.77 0.91 0.60],'LineWidth',2,'DisplayName','loss stable');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');        
           %%
            ax=gca;
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000'};
            %xlim([500 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel('Pupil response - outcome received vs not');
%             set( gca                       , ...
%                 'FontName'   , 'Helvetica' );
%             set([ hXLabel, hYLabel], ...
%                 'FontName'   , 'Helvetica');
% 
%             set([hXLabel, hYLabel]  , ...
%                 'FontSize'   , 13          );

            set(gca, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.02 .02] , ...
                'XMinorTick'  , 'on'      , ...
                'YMinorTick'  , 'on'      , ...
                'YGrid'       , 'off'      , ...
                'XColor'      , [0 0 0 ], ...
                'YColor'      , [0 0 0], ...
                'LineWidth'   , 2         );
            %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8/1.5 6/1.5])
            axis square
            hold on
            legend
            title([group{1},': visit2 vs visit1'])
            saveas(f3,[figdir,'pupil/',group{1},'_win_diff_pupil_response_vol_vs_sta_to_othervol.png']);
            %%
            f4=figure;
            %% plot loss pupil response for each block 
            %  loss bothv-winv
            upper=group_means_diff_delta_loss_bothvvswinv.(group{1})(1,501:end)+group_stes_diff_delta_loss_bothvvswinv.(group{1})(1,501:end);
            lower=group_means_diff_delta_loss_bothvvswinv.(group{1})(1,501:end)-group_stes_diff_delta_loss_bothvvswinv.(group{1})(1,501:end);
            jbfill(xpoints,upper,lower,[0.31 0.15 0.27],[0.31 0.15 0.27],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_delta_loss_bothvvswinv.(group{1})(1,501:end),'Color',[0.31 0.15 0.27],'LineWidth',2,'DisplayName','win volatile');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
            %loss lossv-boths
            upper=group_means_diff_delta_loss_lossvvsboths.(group{1})(1,501:end)+group_stes_diff_delta_loss_lossvvsboths.(group{1})(1,501:end);
            lower=group_means_diff_delta_loss_lossvvsboths.(group{1})(1,501:end)-group_stes_diff_delta_loss_lossvvsboths.(group{1})(1,501:end);
            jbfill(xpoints,upper,lower,[1 0.6 0.60],[1 0.6 0.60],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_delta_loss_lossvvsboths.(group{1})(1,501:end),'Color',[1 0.6 0.60],'LineWidth',2,'DisplayName','win stable');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');        
           %%
            ax=gca;
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000'};
            %xlim([500 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel('Pupil response - outcome received vs not');
%             set( gca                       , ...
%                 'FontName'   , 'Helvetica' );
%             set([ hXLabel, hYLabel], ...
%                 'FontName'   , 'Helvetica');
% 
%             set([hXLabel, hYLabel]  , ...
%                 'FontSize'   , 13          );

            set(gca, ...
                'Box'         , 'off'     , ...
                'TickDir'     , 'out'     , ...
                'TickLength'  , [.02 .02] , ...
                'XMinorTick'  , 'on'      , ...
                'YMinorTick'  , 'on'      , ...
                'YGrid'       , 'off'      , ...
                'XColor'      , [0 0 0 ], ...
                'YColor'      , [0 0 0], ...
                'LineWidth'   , 2         );
            %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8/1.5 6/1.5])
            axis square
            hold on
            legend
            title([group{1},': visit2 vs visit1'])
            saveas(f4,[figdir,'pupil/',group{1},'_loss_diff_pupil_response_vol_vs_sta_to_othervol.png']);
  end
