function [PDC, Direction_Stats,fvec,tsec] = EstimatePDC_STOKV1_Group(Projfolder,varargin)

%% set default values
    opt = ParseArgs(varargin, ...
        'IDs'           ,[],...
        'plotfig'       ,true, ...
        'ModOrds'       ,15,...
        'Cnd'           , 6, ...
        'recalc'        ,false, ...
        'Normalize'     ,'All',...
        'NormPrestim'   ,false,...
        'Freqband'      ,[5 150],...
        'TimeWin'       ,[-200 300],...
        'figpath'       ,fullfile(fileparts(fileparts(Projfolder)),'Results') ...
        );

    animals         =       dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals         =       {animals.name};
    animID          =       cellfun(@(x) x(1:6),animals,'uni',false); % animal IDs
    
    if ~isempty(opt.IDs) % if a subset of animals are selected
        animIdx     =       strcmp(animID,opt.IDs);
        animID      =       animID(animIdx);
        animals     =       animals(animIdx);
    end
    animID{end+1} =     'Average';
    
    opt.Freqband    =       round(opt.Freqband); % Round the frequency band for now
    SaveFileName    =       ['PDCSTOK_Cnd' num2str(opt.Cnd) '_MORD' num2str(opt.ModOrds)];
    ROIs            =       {'V1'};
%% Estimate MVAR parameters using STOK algorithm
    
    if ~exist(fullfile(Projfolder,[SaveFileName '.mat']),'file') || opt.recalc
        for subj    =       1:numel(animals)
            
            %---------------------- Load required Data --------------------
            disp(animID{subj});
            LFP             =       load(fullfile(Projfolder,animals{subj}),'lfpMouse');
            load(fullfile(Projfolder,animals{subj}),'tsec');
            dims            =       size(LFP.lfpMouse);
            LFP.lfpMouse    =       reshape(LFP.lfpMouse, dims(1),dims(2),dims(3)*dims(4),dims(5)); %treat all exps (recordings) as one
            Cond            =       opt.Cnd;
            % (1) Prepare the data
            epochs{subj} = permute(LFP.lfpMouse(:,:,:,Cond), [3,2,1]); %trials, nodes, time
            
            
        end
        % ------------------------ prepare the parameters------------------
        ff          =       .99;
        keepdiag    =       1; % 
        measure     =       'sPDC';
        flow        =       2; % 1 col, 2 row-wise normalization
        fvec        =       round(opt.Freqband(1)):round(opt.Freqband(2));
        load(fullfile(Projfolder,animals{subj}),'srate');
        tvec        =       (tsec>=opt.TimeWin(1)) & (tsec<=opt.TimeWin(2));
        tsec        =       tsec(tvec);
        labels      =       arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
        %----------------- Estimate MVAR params and PDC -------------------
        epochs      =       cellfun(@(x) x(:,:,tvec),epochs,'uni',false); % select the time window 
        % estimate PDC, directionalities using bootstraping
        [pdc,Direction_Stats]  =       bootstrap_PDC(epochs,srate,fvec,tsec,labels,ROIs,...
                            'nboot',        200,...
                            'ModOrds',      opt.ModOrds,...
                            'ff',           ff,...
                            'measure',      measure,...
                            'keepdiag',     keepdiag,...
                            'flow',         flow);
        
        % just prepare PDC structure for the later analysis
        for subj    =       1:numel(animID)
            PDC.(animID{subj})  =       pdc.PDC(:,:,:,:,subj);
            C                   =       pdc.C;
        end


        save(fullfile(Projfolder,SaveFileName),'PDC','C', 'fvec','tsec','labels','Direction_Stats');
    else
        load(fullfile(Projfolder,SaveFileName));
    end
    %% Normalize PDC values using prestim PDCs
    if opt.NormPrestim
        IndNorm     =   find(tsec<0);
        IndNorm     =   IndNorm(100:end); % to remove the unstable part of the PDCs
         for subj   =   1:numel(animals)
            PDC.(animID{subj})  =   PDC.(animID{subj})-mean(PDC.(animID{subj})(:,:,:,IndNorm),4);
            
         end
         
         SaveFigName    =       [SaveFileName '_NormPrestim'];
         ranges         =       [-1 1];   
    else
        SaveFigName     =       SaveFileName;
        ranges          =       [0 1];
    end
    %% Plot the results: individuals and average layer connectivities
    
    if opt.plotfig
        %LFPanalysis.plot_PDC_LFP(PDC,C,tsec,fvec,[-100 300],animID,labels,ranges,ROIs,opt,SaveFigName)
        LFPanalysis.plot_PDC_Directionality(Direction_Stats,tsec,fvec,[-100 300],labels,ROIs,opt,SaveFigName)
    end

end



