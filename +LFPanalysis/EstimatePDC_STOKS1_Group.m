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
    animID{end+1} =     'Average';
    
    opt.Freqband =      round(opt.Freqband); % Round the frequency band for now
    SaveFileName =      ['PDCSTOK_MORD' num2str(opt.ModOrds)];
    ROIs         =      {'cS1','iS1'};
    
%% Estimate MVAR and PDC parameters using STOK algorithm -> 
    % Then compute directionality and apply bootstrap to get the
    % significant effects. All in "bootstrap_PDC" function
    
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
        
        % ------------------------ prepare the parameters------------------
        ff      =       .99;
        keepdiag =      1; % 
        measure =       'sPDC';
        flow    =       2; % 1 col, 2 row-wise normalization
        fvec    =       opt.Freqband(1):opt.Freqband(2);
        load(fullfile(Projfolder,animals{subj}),'srate');
        tvec    =       (tsec>=opt.TimeWin(1)) & (tsec<=opt.TimeWin(2));
        tsec    =       tsec(tvec);
        labels  =       arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
        
        %----------------- Estimate MVAR params and PDC -------------------
        epochs  =       cellfun(@(x) x(:,:,tvec),epochs,'uni',false); % select the time window 
        % estimate PDC, directionalities using bootstraping
        [pdc,Direction_Stats]  =       bootstrap_PDC(epochs,srate,fvec,tsec,...
                            'nboot',        200,...
                            'ModOrds',      opt.ModOrds,...
                            'ff',           ff,...
                            'measure',      measure,...
                            'keepdiag',     keepdiag,...
                            'flow',         flow);
        
        % just prepare PDC structure for the later analysis
        for subj = 1:numel(animID)
            PDC.(animID{subj})  =       pdc(:,:,:,:,subj);
        end
        %--------------------------- Save Results -------------------------
        save(fullfile(Projfolder,SaveFileName),'PDC', 'fvec','tsec','labels','Direction_Stats');
    else
        load(fullfile(Projfolder,SaveFileName));
    end
    
    %% Normalize PDC values using prestim PDCs : just for the sake of visualization -> not included in directionality analysis
    
    if opt.NormPrestim
        IndNorm     =       find(tsec<0);
        IndNorm     =       IndNorm(100:end); % to remove the unstable part of the PDCs
         for subj   =       1:numel(animID)+1
            PDC.(animID{subj})  =       (PDC.(animID{subj})-mean(PDC.(animID{subj})(:,:,:,IndNorm),4));
         end
         
         SaveFigName =      [SaveFileName '_NormPrestim']; % Figure Names
         ranges      =      [-1 1];                       % range for plotting using dynet_connplot
    else
        SaveFigName  =      SaveFileName;
        ranges       =      [0 1];
    end

    %% Plot the results: individuals and average layer connectivities
    
    if opt.plotfig
        
        if ~exist(opt.figpath,'dir') % prepare the folder for the results
            mkdir(opt.figpath);
        end
        
        TimeInd     =       (tsec>=-50) & (tsec<=60); % time window for plotting
        
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
            if true
                for roi = 1:2
                    FIG  =      dynet_connplot(Data((1:6)+(roi-1)*6,(1:6)+(roi-1)*6,:,:),tsec(TimeInd),fvec,labels,ranges, [], [],1);
                    if opt.NormPrestim, colormap(jmaColors('coolhotcortex'));end
                    set(FIG,'unit','inch','position',[0 0 25 20],'color','w')
                    export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_' animID{subj} ROIs{roi}]),'-pdf');
                    close;
                end

            end
            close all;

        end
   
    end
    
    %% Estimate UPWARD and DOWNWARD FCs
    
    Data        =       PDC.Average;
    for ch      =       1:size(Data,1)
        Data(ch,ch,:,:) = 0;
    end
    %% (1) node-wise analysis
    load('LayerColors.mat');
    LNames = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
    %-----------------------Nodes outflows-----------------------------------
    if opt.plotfig
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
                xlim([-50 60])
                if opt.NormPrestim
                    %ylim([-0.4 0.5])
                else
                    %ylim([0 .7])
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

    % -----------------------Average Outflows---------------------------------
    if opt.plotfig
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
            xlim([-50 60])
            xlabel('Time (ms)')
            vline([0],'k--')
            legend(SP,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));
            set(gca,'fontsize',16);
            set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
            title ('Average Outflow of layers')
            export_fig(Fig2,fullfile(opt.figpath,[SaveFigName '_AveragewOutflow_' ROIs{roi}]),'-pdf');
        end
    end
    close all;
    
    
%% ---------------------------Directionality-----------------------------
FS  =   12;
if  opt.plotfig   
    for roi     =   1:2 % for each hemi
       lnum     =   size(Data,1)/2;
              
       %-------------------------Averaged Over Frequencies-----------------
        Fig2 = figure;
        clear SP;
        line([-50 100],[0 0],'linestyle','--','color','k','linewidth',1.3);
        hold on;
        for i = 2:lnum-1
            LIdx = find(strcmp({Direction_Stats.LName},[labels{i} '_' ROIs{roi}]));
            SP(i-1) = plot(tsec,sum(Direction_Stats(LIdx).Median,1),'color',Colors(i,:),'linewidth',2);
        end
        xlim([-50 60])
        xlabel('Time (ms)')
        vline([0],'k--')
        legend(SP,labels(2:lnum-1));
        set(gca,'fontsize',16);
        set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
        title ('Upward flow of layers')
        export_fig(Fig2,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AveragedFreqs_' ROIs{roi}]),'-pdf');

        %--------------------------All Frequencies-------------------------
        ind         =       (tsec>=-50) & (tsec<=60); % time window for plotting
        tseccorr    =       tsec(ind);
        Fig22       =       figure;
        tsec2       =       (tsec>=0);
        for i = 2:lnum-1
            
            SP(i-1)     =       subplot(4,1,i-1);
            LIdx        =       find(strcmp({Direction_Stats.LName},[labels{i} '_' ROIs{roi}]));
            Dataplot    =       Direction_Stats(LIdx).Median(:,ind);
            Transplot   =       Direction_Stats(LIdx).PostStim.h;
            TP = zeros(size(Direction_Stats(LIdx).Median));
            TP(:,tsec2) = Transplot;
            TP = TP(:,ind);
            imagesc(Dataplot,'alphadata',TP*.5+.5);
            axis xy;
            xtickso = find(ismember(tseccorr,-100:20:100));
            set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',1:40:numel(fvec),'yticklabel',fvec(1:40:end));
            M = max(abs(Dataplot(:)));
            caxis([-M M]);
            %caxis([-.2 .2]);
            vline(xtickso(tseccorr(xtickso)==0),'k--');
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
        
        LIdx        =       find(strcmp({Direction_Stats.LName},['All_' ROIs{roi}]));
        Dataplot    =       Direction_Stats(LIdx).Median(:,ind);
        Transplot   =       Direction_Stats(LIdx).PostStim.h;
        TP = zeros(size(Direction_Stats(LIdx).Median));
        TP(:,tsec2) = Transplot;
        TP = TP(:,ind);
        
        Fig3 = figure;
        imagesc(Dataplot,'alphadata',TP*.5+.5);
        axis xy;
        xtickso = find(ismember(tseccorr,-100:20:100));
        set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',1:40:numel(fvec),'yticklabel',fvec(1:40:end));
        M = max(abs(Dataplot(:)));
        caxis([-M M]);
        vline(xtickso(tseccorr(xtickso)==0),'k--');
        colorbar;
        xlabel('Time(msec)');
        ylabel('Frequency(Hz)');
        set(gca,'fontsize',FS);
        colormap(jmaColors('coolhotcortex'))
        set(Fig3,'unit','inch','position',[1 1 10 4],'color','w');
        export_fig(Fig3,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AveragedLayers_' ROIs{roi}]),'-pdf');

        close all

    end

end

