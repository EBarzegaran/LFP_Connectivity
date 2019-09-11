function [PDC, fvec] = EstimatePDC_STOKS1_Group(Projfolder,varargin)


%% set default values
    opt = ParseArgs(varargin, ...
        'plotfig'       ,true, ...
        'ModOrds'       ,15,...
        'Cnd'           ,6, ...
        'recalc'        ,false, ...
        'Normalize'     ,'All',...
        'Freqband'      ,[20 150],...
        'SStat'         ,'Paired',...
        'doStats'       ,false,...
        'NormPrestim'   ,false,...
        'doPermute'     ,false,...
        'typePermute'   ,'layers',...
        'NPerm'         ,100,...
        'TimeWin'       ,[-100 100],...
        'figpath'       ,fullfile(Projfolder,'Results') ...
        );

    animals     =       dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals     =       {animals.name};
    animID      =       cellfun(@(x) x(1:4),animals,'uni',false); % animal IDs
    
    opt.Freqband =      round(opt.Freqband); % Round the frequency band for now
    SaveFileName =      ['PDCSTOK_MORD' num2str(opt.ModOrds)];
    ROIs         =      {'cS1','iS1'};
    
%% Estimate MVAR and PDC parameters using STOK algorithm
    
    if ~exist(fullfile(Projfolder,[SaveFileName '.mat']),'file') || opt.recalc
        disp('Reading LFP data');
        for subj = 1:numel(animals)
            
            %---------------------- Load the Data -------------------------
            disp(animID{subj});
            LFP     =       load(fullfile(Projfolder,animals{subj}),'lfpRat');
            load(fullfile(Projfolder,animals{subj}),'tsec');
           
            % --------------------- Prepare the data ----------------------
            epochs{subj}    =       permute(LFP.lfpRat(:,:,:), [3,2,1]); %trials, nodes, time
        end
        
        %Epochs  =       cat(1,epochs{:}); % Agreggate data for PDC estimation

        % ------------------------ prepare the parameters------------------
        ff      =       .99;
        keepdiag =      1; % 
        measure =       'PDCnn';
        flow    =       2; % 1 col, 2 row-wise normalization
        fvec    =       opt.Freqband(1):opt.Freqband(2);
        load(fullfile(Projfolder,animals{subj}),'srate');
        tvec    =       (tsec>=opt.TimeWin(1)) & (tsec<=opt.TimeWin(2));
        tsec    =       tsec(tvec);
        labels  =       arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
        
        %----------------- Estimate MVAR params and PDC -------------------
        epochs  =       cellfun(@(x) x(:,:,tvec),epochs,'uni',false); % select the time window 
        [pdc,Direction_Stats] =       bootstrap_PDC(epochs,srate,fvec,tsec,...
                            'nboot',        200,...
                            'ModOrds',      opt.ModOrds,...
                            'ff',           ff,...
                            'measure',      measure,...
                            'keepdiag',     keepdiag,...
                            'flow',         flow);
        
        % just prepare PDC structure for the later analysis
        PDC.(animID{S})  =       pdc(:,:,:,:,S);
        %--------------------------- Save Results -------------------------
        save(fullfile(Projfolder,SaveFileName),'PDC', 'fvec','tsec','labels','Direction_Stats');
    else
        load(fullfile(Projfolder,SaveFileName));
    end
    
    %% Normalize PDC values using prestim PDCs
    
    if opt.NormPrestim
        IndNorm     =       find(tsec<0);
        IndNorm     =       IndNorm(100:end); % to remove the unstable part of the PDCs
         for subj   =       1:numel(animals)
            PDC.(animID{subj})  =       (PDC.(animID{subj})-mean(PDC.(animID{subj})(:,:,:,IndNorm),4));
            
         end
         
         SaveFigName =      [SaveFileName '_NormPrestim'];
    else
        SaveFigName  =      SaveFileName;
    end

    %% Plot the results: individuals and average layer connectivities
    
    
    if opt.plotfig
        if ~exist(opt.figpath,'dir')
            mkdir(opt.figpath);
        end
        
        for subj    =       1:numel(animals)
            Data    =       PDC.(animID{subj})(:,:,:,:);
            % what kind of normalization to be done
            if strcmpi(opt.Normalize,'Channel')% (1) normalize over every single channel
                Data =      Data./repmat(max(max(Data,[],4),[],3),[1 1 size(Data,3) size(Data,4)]); 
                
            elseif strcmpi(opt.Normalize,'All')% (2) normalize over all channels
                for ch =    1:size(Data,1)
                    Data(ch,ch,:,:) = 0;
                end
                Data =      Data./max(Data(:));
            end
            %plot individual animal resaults
            if false
                FIG  =      dynet_connplot(Data(1:6,1:6,:,:),tsec,fvec,labels, [], [], [],0);
                set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
                export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_' animID{subj} '_cS1']),'-pdf');
                close;

                FIG  =      dynet_connplot(Data(7:12,7:12,:,:),tsec,fvec,labels, [], [], [],0);
                set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
                export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_' animID{subj} '_cS1']),'-pdf');
                close;
            end
            close all;
            
            DataM(:,:,:,:,subj)     =   Data;%(Data-Databm)./Databs;
        end
        
        % plot average over all animals
        if opt.NormPrestim
            ranges      =       [-1 1];
        else
            ranges      =       [0 1];
        end
        TimeInd         =       (tsec>=-50);
        Data            =       (mean(DataM(1:6,1:6,:,TimeInd,:),5));
        FIG             =       dynet_connplot(Data./max(Data(:)),    tsec(TimeInd),  fvec,   labels, ranges, [], [],0);
        if opt.NormPrestim, colormap(jmaColors('coolhotcortex'));end
        set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
        export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_ AverageAll_' ROIs{1}]),'-pdf');
        
        Data            =       mean(DataM(7:12,7:12,:,TimeInd,:),5);
        FIG             =       dynet_connplot(Data./max(Data(:)),    tsec(TimeInd), fvec,    labels, ranges, [], [],0);
        if opt.NormPrestim, colormap(jmaColors('coolhotcortex'));end
        set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
        export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_ AverageAll_' ROIs{2}]),'-pdf');
        
              
    end
    
    %% estimate UPWARD and DOWNWARD FCs
    clear DataM
    for subj    =      1:numel(animals)% Aggregate animals data        
        Data    =      PDC.(animID{subj})(:,:,:,:);
        for ch  =      1:size(Data,1)
            Data(ch,ch,:,:) = 0;
        end
        Data    =      Data./max(Data(:));      
        DataM(:,:,:,:,subj) =   Data;
    end 
    
    Data        =       mean(DataM(:,:,:,:,:),5);
   
    %% (1) node-wise analysis
    load('LayerColors.mat');
    LNames = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
    %-----------------------Nodes outflows-----------------------------------
    if true
        lnum  = size(Data,1)/2;
        for roi = 1:2
            Fig1 = figure;
            for i = 1:lnum
                subplot(lnum,1,i);
                for j  = 1:lnum
                    SP(j)=plot(tsec,squeeze(mean(Data(j+(roi-1)*lnum,i+(roi-1)*lnum,:,:),3)),'color',Colors(j,:),'linewidth',2);
                    hold on;
                end
                title(['L' num2str(i)]);
                xlim([-50 100])
                if opt.NormPrestim
                    ylim([-0.4 0.5])
                else
                    ylim([0 .7])
                end
                vline([0],'k--')
            end
            legend(SP,LNames);
            xlabel('Time (ms)');
            ylabel('outPDC')
            set(Fig1,'unit','inch','position',[2 0 15 20],'color','w')
            export_fig(Fig1,fullfile(opt.figpath,[SaveFigName '_LayersOutflow_' ROIs{roi}]),'-pdf');
        end
    end

    %% -----------------------Average Outflows---------------------------------
    if true
        for roi = 1:2
            lnum  = size(Data,1)/2;
            Data2 = Data((1:lnum)+(roi-1)*lnum,(1:lnum)+(roi-1)*lnum,:,:);
            PDCavgOut = arrayfun(@(x) squeeze(mean(Data2([1:x-1 x+1:end],x,:,:),1)),1:size(Data2,1),'uni',false);

            Fig2 = figure;
            line([-100 300],[0 0],'linestyle','--','color','k','linewidth',1.3);
            hold on;
            for i = 1:numel(PDCavgOut)
                SP(i) = plot(tsec,mean(PDCavgOut{i},1),'color',Colors(i,:),'linewidth',2);
            end
            xlim([-50 100])
            xlabel('Time (ms)')
            vline([0],'k--')
            legend(SP,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));
            set(gca,'fontsize',16);
            set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
            title ('Average Outflow of layers')
            export_fig(Fig2,fullfile(opt.figpath,[SaveFigName '_AveragewOutflow_' ROIs{roi}]),'-pdf');
        end
    end
    %close all;
    %% ---------------------------Directionality-----------------------------
    FS = 12;
    
    for roi =1:2
       lnum  = size(DataM,1)/2;
       DataM2 = abs(DataM((1:lnum)+(roi-1)*lnum,(1:lnum)+(roi-1)*lnum,:,:,:)); 
       %Diff = arrayfun(@(x) squeeze(mean(DataM2(1:x-1,x,:,:,:),1) - mean(DataM2(x+1:end,x,:,:,:),1))./squeeze(mean(cat(1,mean(DataM2(1:x-1,x,:,:,:),1),mean(DataM2(x+1:end,x,:,:,:),1)),1)),1:6,'uni',false);
       Diff = arrayfun(@(x) squeeze(mean(DataM2(1:x-1,x,:,:,:),1) - mean(DataM2(x+1:end,x,:,:,:),1)),1:6,'uni',false);%./squeeze(mean(Data(x,x:end,:,:),2) + mean(Data(x,1:x,:,:),2)),1:size(Data,1),'uni',false);
       Diff = cellfun(@(x) mean(x,3),Diff,'uni',false);
       
       if true
           %-------------------------Averaged Over Frequencies-----------------
            Fig2 = figure;
            line([-50 100],[0 0],'linestyle','--','color','k','linewidth',1.3);
            hold on;
            for i = 2:numel(PDCavgOut)-1
                SP(i-1) = plot(tsec,sum(Diff{i},1),'color',Colors(i,:),'linewidth',2);
            end
            xlim([-50 100])
            xlabel('Time (ms)')
            vline([0],'k--')
            legend(SP,arrayfun(@(x) ['L' num2str(x)],2:5,'uni',false));
            set(gca,'fontsize',16);
            set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
            title ('Upward flow of layers')
            export_fig(Fig2,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AveragedFreqs_' ROIs{roi}]),'-pdf');

            %--------------------------All Frequencies-------------------------
            ind = (tsec>=-100);
            tseccorr = tsec(ind);
            Fig22 = figure;
            for i = 2:numel(PDCavgOut)-1
                SP(i-1) = subplot(4,1,i-1);
                Dataplot = Diff{i}(:,ind);
                imagesc(Dataplot);
                axis xy;
                xtickso = find(ismember(tseccorr,-100:100:300));
                set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',1:40:numel(fvec),'yticklabel',fvec(1:40:end));
                M = max(abs(Diff{i}(:)));
                %caxis([-M M]);
                caxis([-.35 .35]);
                vline(xtickso(2),'k--');
                colorbar;
                if i==5
                    xlabel('Time(msec)');
                    ylabel('Frequency(Hz)');
                end
                title(labels{i})
                set(gca,'fontsize',FS)
            end
            colormap(jmaColors('coolhotcortex'))
            set(Fig22,'unit','inch','position',[1 1 10 15],'color','w')

            export_fig(Fig22,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AllFreqs_' ROIs{roi}]),'-pdf');

            %--------------------average over all Layers-----------------------
            DiffM = mean(cat(3,Diff{2:5}),3);
            Fig3 = figure;
            imagesc(DiffM(:,ind));
            axis xy;
            xtickso = find(ismember(tseccorr,-100:100:300));
            set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',1:40:numel(fvec),'yticklabel',fvec(1:40:end));
            caxis([-.2 .2]);
            vline(xtickso(2),'k--');
            colorbar;
            xlabel('Time(msec)');
            ylabel('Frequency(Hz)');
            set(gca,'fontsize',FS);
            colormap(jmaColors('coolhotcortex'))
            set(Fig3,'unit','inch','position',[1 1 10 4],'color','w');
            export_fig(Fig3,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AveragedLayers_' ROIs{roi}]),'-pdf');

            %close all
       end

    end

end

