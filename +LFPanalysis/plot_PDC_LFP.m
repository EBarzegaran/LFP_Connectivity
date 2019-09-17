function plot_PDC_LFP(PDC,C,tsec,fvec,tWin,animID,labels,ranges,ROIs,opt,SaveFigName)


% This function needs to be cleaned
%% PDC results
        if ~exist(opt.figpath,'dir') % prepare the folder for the results
            mkdir(opt.figpath);
        end
        
        TimeInd     =       (tsec>=tWin(1)) & (tsec<=tWin(2)); % time window for plotting
        
        for subj    =       numel(animID):numel(animID)
            Data    =       PDC.(animID{subj})(:,:,:,TimeInd);
            
            %---------- what kind of normalization to be done -------------
            if strcmpi(opt.Normalize,'Channel')% (1) normalize over every single channel
                Data    =      Data./repmat(max(max(Data,[],4),[],3),[1 1 size(Data,3) size(Data,4)]); 
                
            elseif strcmpi(opt.Normalize,'All')% (2) normalize over all channels
                for ch  =      1:size(Data,1)
                    Data(ch,ch,:,:) = 0;
                end
                Data    =      Data./max(Data(:));
            end
            
            %------------plot individual animal resaults-------------------
            for roi = 1:numel(ROIs)
                FIG  =      dynet_connplot(Data((1:6)+(roi-1)*6,(1:6)+(roi-1)*6,:,:),tsec(TimeInd),fvec,labels,ranges*.996, [], [],1);
                if opt.NormPrestim, colormap(jmaColors('coolhotcortex'));end
                set(FIG,'unit','inch','position',[0 0 25 20],'color','w')
                export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_' animID{subj} ROIs{roi}]),'-pdf');
                close;
            end
        end
            %-------------- plot C tuning parameter -----------------------
            FIG = figure;
            subplot(2,1,1); 
            hold on;box on;
            load LayerColors.mat;
            for l = 1:size(C,2)
                plot(tsec,C(:,l),'color',Colors(l,:),'linewidth',1.5)
            end
            vline(0,'k--')
            xlim([tWin(1) tWin(2)])
            legend(animID)
            title('C Parameter')
            set(gca,'fontsize',12)
            
            subplot(2,1,2);
            plot(tsec,mean(C,2),'color','k','linewidth',1.5);
            xlabel('Time(mS)')
            title('Averaged C')
            vline(0,'k--')
            xlim([tWin(1) tWin(2)])
            set(FIG,'unit','inch','position',[0 0 15 10],'color','w')
            set(gca,'fontsize',12)
            export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_C_Param']),'-pdf');
            
            close all;
        
    Data        =       PDC.Average;
    for ch      =       1:size(Data,1)
        Data(ch,ch,:,:) = 0;
    end
%% (1) node-wise outflow results
    load('LayerColors.mat');
    LNames = labels;
    %-----------------------Nodes outflows-----------------------------------
    lnum  = size(Data,1)/numel(ROIs);
    for roi = 1:numel(ROIs)
        Fig1 = figure;
        for i = 1:lnum
            subplot(lnum,1,i);
            for j  = 1:lnum
                SP(j)=plot(tsec,squeeze(mean(Data(j+(roi-1)*lnum,i+(roi-1)*lnum,:,:),3)),'color',Colors(j,:),'linewidth',2);
                hold on;
            end
            title(['L' num2str(i)]);
            xlim([tWin(1) tWin(2)])
            vline([0],'k--')
        end
        legend(SP,LNames);
        xlabel('Time (ms)');
        ylabel('outPDC')
        set(Fig1,'unit','inch','position',[2 0 15 20],'color','w')
        export_fig(Fig1,fullfile(opt.figpath,[SaveFigName '_LayersOutflow_' ROIs{roi}]),'-pdf');
    end

    % -----------------------Average Outflows---------------------------------
    for roi = 1:numel(ROIs)
        Data2       =       Data((1:lnum)+(roi-1)*lnum,(1:lnum)+(roi-1)*lnum,:,:);
        PDCavgOut   =       arrayfun(@(x) squeeze(mean(Data2([1:x-1 x+1:end],x,:,:),1)),1:size(Data2,1),'uni',false);

        Fig2        =       figure;
        line(tWin,[0 0],'linestyle','--','color','k','linewidth',1.3);
        hold on;
        for i       =       1:numel(PDCavgOut)
            SP(i)   =       plot(tsec,mean(PDCavgOut{i},1),'color',Colors(i,:),'linewidth',2);
        end
        xlim([tWin(1) tWin(2)])
        xlabel('Time (ms)')
        vline([0],'k--')
        legend(SP,arrayfun(@(x) ['L' num2str(x)],1:lnum,'uni',false));
        set(gca,'fontsize',16);
        set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
        title ('Average Outflow of layers')
        export_fig(Fig2,fullfile(opt.figpath,[SaveFigName '_AveragewOutflow_' ROIs{roi}]),'-pdf');
    end
    close all;
end