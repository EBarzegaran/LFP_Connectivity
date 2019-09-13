function [PDC, fvec] = EstimatePDC_STOKS1_Group(Projfolder,varargin)


%% set default values
    opt = ParseArgs(varargin, ...
        'plotfig'       ,true, ...
        'ModOrds'       ,15,...
        'Cnd'           ,6, ...
        'recalc'        ,false, ...
        'StatSide'      ,'right',...
        'Normalize'     ,'All',...
        'Freqband'      ,[20 150],...
        'NormPrestim'   ,false,...
        'TimeWin'       ,[-100 100],...
        'figpath'       ,fullfile(Projfolder,'Results') ...
        );

    animals     =       dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals     =       {animals.name};
    animID      =       cellfun(@(x) x(1:4),animals,'uni',false); % animal IDs
    animID{end+1} =     'Average';
    
    opt.Freqband =      round(opt.Freqband); % Round the frequency band for now
    SaveFileName =      ['PDCSTOK_MORD' num2str(opt.ModOrds) '_' opt.StatSide];
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
        [pdc,Direction_Stats]  =       bootstrap_PDC(epochs,srate,fvec,tsec, labels, ROIs,...
                            'nboot',        200,...
                            'ModOrds',      opt.ModOrds,...
                            'ff',           ff,...
                            'measure',      measure,...
                            'keepdiag',     keepdiag,...
                            'flow',         flow,...
                            'StatSide',     opt.StatSide);
        
        % just prepare PDC structure for the later analysis
        for subj =      1:numel(animID)
            PDC.(animID{subj})  =       pdc.PDC(:,:,:,:,subj);
            C                   =       pdc.C;
        end
        %--------------------------- Save Results -------------------------
        save(fullfile(Projfolder,SaveFileName),'PDC','C', 'fvec','tsec','labels','Direction_Stats');
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
        LFPanalysis.plot_PDC_LFP(PDC,C,tsec,fvec,[-50 60],animID,labels,ranges,ROIs,opt,SaveFigName)
        LFPanalysis.plot_PDC_Directionality(Direction_Stats,tsec,fvec,[-50 60],labels,ROIs,opt,SaveFigName)
    end

end

