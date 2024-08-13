%% plot winvol winstable lossvol lossstable across blocks\
for group={'reboxetine','placebo'}
      for visit=1:2

            f1=figure;%main figure combining both time courses [reaction to chosen rewards]
            xpoints=1:3000;
            y(1,1:length(xpoints))=zeros; %%reference line
           
            %% plot win pupil response for each block 
            %  both vol block
            upper=reshape(group_means_rew_delivered.(group{1})(visit,1,501:end)+group_stes_rew_delivered.(group{1})(visit,1,501:end),1,3000);
            lower=reshape(group_means_rew_delivered.(group{1})(visit,1,501:end)-group_stes_rew_delivered.(group{1})(visit,1,501:end),1,3000);
            jbfill(xpoints,upper,lower,[.15 0.35 0.26],[0 1 0],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_rew_delivered.(group{1})(visit,1,501:end),1,3000),'Color',[.15 0.35 0.26],'LineWidth',2,'DisplayName','win both vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
            %win vol block
            upper=reshape(group_means_rew_delivered.(group{1})(visit,2,501:end)+group_stes_rew_delivered.(group{1})(visit,2,501:end),1,3000);
            lower=reshape(group_means_rew_delivered.(group{1})(visit,2,501:end)-group_stes_rew_delivered.(group{1})(visit,2,501:end),1,3000);
            jbfill(xpoints,upper,lower,[0.46 .68 0.61],[.45 0.5 0.55],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_rew_delivered.(group{1})(visit,2,501:end),1,3000),'Color',[0.46 .68 0.61],'LineWidth',2,'DisplayName','win win vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            %loss vol block
            upper=reshape(group_means_rew_delivered.(group{1})(visit,3,501:end)+group_stes_rew_delivered.(group{1})(visit,3,501:end),1,3000);
            lower=reshape(group_means_rew_delivered.(group{1})(visit,3,501:end)-group_stes_rew_delivered.(group{1})(visit,3,501:end),1,3000);
            jbfill(xpoints,upper,lower,[0.55 .85 0.55],[.45 0.5 0.55],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_rew_delivered.(group{1})(visit,3,501:end),1,3000),'Color',[0.55 .85 0.55],'LineWidth',2,'DisplayName','win loss vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            %both stable block
            upper=reshape(group_means_rew_delivered.(group{1})(visit,4,501:end)+group_stes_rew_delivered.(group{1})(visit,4,501:end),1,3000);
            lower=reshape(group_means_rew_delivered.(group{1})(visit,4,501:end)-group_stes_rew_delivered.(group{1})(visit,4,501:end),1,3000);
            jbfill(xpoints,upper,lower,[0.77 0.91 0.60],[.45 0.5 0.55],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_rew_delivered.(group{1})(visit,4,501:end),1,3000),'Color',[0.77 0.91 0.60],'LineWidth',2,'DisplayName','win both stable');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');        
           %%
            ax=gca;
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000'};
            %xlim([500 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel('Pupil response - outcome received');
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
            saveas(f1,[figdir,'pupil/',group{1},'_visit',num2str(visit),'_win_delivered_for_each_block.png']);
            %% -------
            f2=figure;
            %xpoints=1:length(lossvol_lossv);
            y(1,1:length(xpoints))=zeros; %%reference line

            %% plot loss pupil response for each block 
            %  both vol block
            upper=reshape(group_means_pun_delivered.(group{1})(visit,1,501:end)+group_stes_pun_delivered.(group{1})(visit,1,501:end),1,3000);
            lower=reshape(group_means_pun_delivered.(group{1})(visit,1,501:end)-group_stes_pun_delivered.(group{1})(visit,1,501:end),1,3000);
            jbfill(xpoints,upper,lower,[0.31 0.15 0.27],[0.31 0.15 0.27],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_pun_delivered.(group{1})(visit,1,501:end),1,3000),'Color',[0.31 0.15 0.27],'LineWidth',2,'DisplayName','loss both vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
            %loss vol block
            upper=reshape(group_means_pun_delivered.(group{1})(visit,2,501:end)+group_stes_pun_delivered.(group{1})(visit,2,501:end),1,3000);
            lower=reshape(group_means_pun_delivered.(group{1})(visit,2,501:end)-group_stes_pun_delivered.(group{1})(visit,2,501:end),1,3000);
            jbfill(xpoints,upper,lower,[0.91 .27 0.33],[0.91 .27 0.33],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_pun_delivered.(group{1})(visit,2,501:end),1,3000),'Color',[0.91 .27 0.33],'LineWidth',2,'DisplayName','loss win vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            %loss vol block
            upper=reshape(group_means_pun_delivered.(group{1})(visit,3,501:end)+group_stes_pun_delivered.(group{1})(visit,3,501:end),1,3000);
            lower=reshape(group_means_pun_delivered.(group{1})(visit,3,501:end)-group_stes_pun_delivered.(group{1})(visit,3,501:end),1,3000);
            jbfill(xpoints,upper,lower,[0.65 0.1 0.18],[0.65 0.1 0.18],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_pun_delivered.(group{1})(visit,3,501:end),1,3000),'Color',[0.65 0.1 0.18],'LineWidth',2,'DisplayName','loss loss vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            %both stable block
            upper=reshape(group_means_pun_delivered.(group{1})(visit,4,501:end)+group_stes_pun_delivered.(group{1})(visit,4,501:end),1,3000);
            lower=reshape(group_means_pun_delivered.(group{1})(visit,4,501:end)-group_stes_pun_delivered.(group{1})(visit,4,501:end),1,3000);
            jbfill(xpoints,upper,lower,[1 0.6 0.60],[1 0.6 0.60],0, 0.15);
            hold on
            plot(xpoints,reshape(group_means_pun_delivered.(group{1})(visit,4,501:end),1,3000),'Color',[1 0.6 0.60],'LineWidth',2,'DisplayName','loss both stable');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');        
          %%
            ax=gca;
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000'};
            %xlim([500 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel('Pupil response - outcome received');
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

            saveas(f2,[figdir,'pupil/',group{1},'_visit',num2str(visit),'_loss_delivered_for_each_block.png']);
      end
    %% plot win vol visit2 - visit1
            f3=figure;
            %% plot win pupil response for each block 
            %  both vol block
            upper=group_means_diff_rew_delivered.(group{1})(1,501:end)+group_stes_diff_rew_delivered.(group{1})(1,501:end);
            lower=group_means_diff_rew_delivered.(group{1})(1,501:end)-group_stes_diff_rew_delivered.(group{1})(1,501:end);
            jbfill(xpoints,upper,lower,[.15 0.35 0.26],[0 1 0],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_rew_delivered.(group{1})(1,501:end),'Color',[.15 0.35 0.26],'LineWidth',2,'DisplayName','win both vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
            %win vol block
            upper=group_means_diff_rew_delivered.(group{1})(2,501:end)+group_stes_diff_rew_delivered.(group{1})(2,501:end);
            lower=group_means_diff_rew_delivered.(group{1})(2,501:end)-group_stes_diff_rew_delivered.(group{1})(2,501:end);
            jbfill(xpoints,upper,lower,[0.46 .68 0.61],[.45 0.5 0.55],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_rew_delivered.(group{1})(2,501:end),'Color',[0.46 .68 0.61],'LineWidth',2,'DisplayName','win win vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            %loss vol block
            upper=group_means_diff_rew_delivered.(group{1})(3,501:end)+group_stes_diff_rew_delivered.(group{1})(3,501:end);
            lower=group_means_diff_rew_delivered.(group{1})(3,501:end)-group_stes_diff_rew_delivered.(group{1})(3,501:end);
            jbfill(xpoints,upper,lower,[0.55 .85 0.55],[.45 0.5 0.55],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_rew_delivered.(group{1})(3,501:end),'Color',[0.55 .85 0.55],'LineWidth',2,'DisplayName','win loss vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            %both stable block
            upper=group_means_diff_rew_delivered.(group{1})(4,501:end)+group_stes_diff_rew_delivered.(group{1})(4,501:end);
            lower=group_means_diff_rew_delivered.(group{1})(4,501:end)-group_stes_diff_rew_delivered.(group{1})(4,501:end);
            jbfill(xpoints,upper,lower,[0.77 0.91 0.60],[.45 0.5 0.55],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_rew_delivered.(group{1})(4,501:end),'Color',[0.77 0.91 0.60],'LineWidth',2,'DisplayName','win both stable');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');        
           %%
            ax=gca;
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000'};
            %xlim([500 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel('pupil response - received');
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
            saveas(f3,[figdir,'pupil/',group{1},'_win_diff_delivered_for_each_block.png']);
            %%
            f4=figure;
            %% plot loss pupil response for each block 
            %  both vol block
            upper=group_means_diff_pun_delivered.(group{1})(1,501:end)+group_stes_diff_pun_delivered.(group{1})(1,501:end);
            lower=group_means_diff_pun_delivered.(group{1})(1,501:end)-group_stes_diff_pun_delivered.(group{1})(1,501:end);
            jbfill(xpoints,upper,lower,[0.31 0.15 0.27],[0.31 0.15 0.27],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_pun_delivered.(group{1})(1,501:end),'Color',[0.31 0.15 0.27],'LineWidth',2,'DisplayName','loss both vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            hold on
            %loss vol block
            upper=group_means_diff_pun_delivered.(group{1})(2,501:end)+group_stes_diff_pun_delivered.(group{1})(2,501:end);
            lower=group_means_diff_pun_delivered.(group{1})(2,501:end)-group_stes_diff_pun_delivered.(group{1})(2,501:end);
            jbfill(xpoints,upper,lower,[0.91 .27 0.33],[0.91 .27 0.33],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_pun_delivered.(group{1})(2,501:end),'Color',[0.91 .27 0.33],'LineWidth',2,'DisplayName','loss win vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            %loss vol block
            upper=group_means_diff_pun_delivered.(group{1})(3,501:end)+group_stes_diff_pun_delivered.(group{1})(3,501:end);
            lower=group_means_diff_pun_delivered.(group{1})(3,501:end)-group_stes_diff_pun_delivered.(group{1})(3,501:end);
            jbfill(xpoints,upper,lower,[0.65 0.1 0.18],[0.65 0.1 0.18],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_pun_delivered.(group{1})(3,501:end),'Color',[0.65 0.1 0.18],'LineWidth',2,'DisplayName','loss loss vol');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');
            %both stable block
            upper=group_means_diff_pun_delivered.(group{1})(4,501:end)+group_stes_diff_pun_delivered.(group{1})(4,501:end);
            lower=group_means_diff_pun_delivered.(group{1})(4,501:end)-group_stes_diff_pun_delivered.(group{1})(4,501:end);
            jbfill(xpoints,upper,lower,[1 0.6 0.60],[1 0.6 0.60],0, 0.15);
            hold on
            plot(xpoints,group_means_diff_pun_delivered.(group{1})(4,501:end),'Color',[1 0.6 0.60],'LineWidth',2,'DisplayName','loss both stable');
            hold on
            plot(xpoints,y,'--k','LineWidth',1,'HandleVisibility','off');        
           %%
            ax=gca;
            ax.XTickLabel = {'-1000','0','1000','2000','3000','4000','5000'};
            %xlim([500 3500])
            set(gca,'xgrid','on')
            hXLabel=xlabel('Time in ms from delivery of outcome');
            hYLabel=ylabel('pupil response - received');
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
            saveas(f4,[figdir,'pupil/',group{1},'_loss_diff_delivered_for_each_block.png']);
  end
