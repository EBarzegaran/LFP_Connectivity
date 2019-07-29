function [PDC, fvec] = EstimatePDC_STOK(Projfolder,varargin)

%% set default values
    opt = ParseArgs(varargin, ...
        'plotfig'  ,true, ...
        'ModOrds'  ,15,...
        'Cnd'      , 6, ...
        'recalc'   ,false, ...
        'Normalize','All',...
        'Freqband' , [5 150],...
        'TimeWin'  , [-300 300],...
        'figpath'  ,fullfile(fileparts(fileparts(Projfolder)),'Results') ...
        );

    animals = dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals = {animals.name};
    animID = cellfun(@(x) x(1:6),animals,'uni',false); % animal IDs
    
    opt.Freqband = round(opt.Freqband); % Round the frequency band for now
    
%% Estimate MVAR parameters using STOK algorithm
    xticklabels = -500:100:500;
    xticks = linspace (0,1250,numel(xticklabels));
    
    if ~exist(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd) '.mat']),'file') || opt.recalc
        for subj = 1:numel(animals)
            
            %---------------------- Load required Data --------------------
            disp(animID{subj});
            LFP = load(fullfile(Projfolder,animals{subj}),'lfpMouse');
            load(fullfile(Projfolder,animals{subj}),'tsec');
            dims = size(LFP.lfpMouse);
            LFP.lfpMouse = reshape(LFP.lfpMouse, dims(1),dims(2),dims(3)*dims(4),dims(5)); %treat all exps (recordings) as one

            %--------------- Estimate MVAR params and PDC -----------------
            Cond = opt.Cnd;
            % (1) Prepare the data
            epochs = permute(LFP.lfpMouse(:,:,:,Cond), [3,2,1]); %trials, nodes, time
            % (2) prepare the parameters
            ff=.99;
            keepdiag = 1; % 
            measure='sPDC';
            flow = 2; % 1 col, 2 row-wise normalization
            fvec=1:150;
            load(fullfile(Projfolder,animals{subj}),'srate');
            tvec = (tsec>=opt.TimeWin(1)) & (tsec<=opt.TimeWin(2));
            tsec = tsec(tvec);
            labels = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
            % (3) estimate MVAR coeffs using stok
            KF = dynet_SSM_STOK(epochs(:,:,tvec),opt.ModOrds,ff);  
            % (4) estimate PDC based on MVAR coeffs
            PDC.(animID{subj}) = dynet_ar2pdc(KF,srate,fvec,measure,keepdiag,flow);
        end

        save(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd)]),'PDC', 'fvec','tsec','labels');
    else
        load(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd)]));
    end

    %% Plot the results: individuals and average layer connecdtivities
    
    if opt.plotfig
        if ~exist(opt.figpath,'dir')
            mkdir(opt.figpath);
        end
        
        for subj = 1:numel(animals)
            Data = PDC.(animID{subj})(:,:,opt.Freqband(1):opt.Freqband(2),:);
            % what kind of normalization to be done
            if strcmpi(opt.Normalize,'Channel')% (1) normalize over every single channel
                Data = Data./repmat(max(max(Data,[],4),[],3),[1 1 size(Data,3) size(Data,4)]); 
                
            elseif strcmpi(opt.Normalize,'All')% (2) normalize over all channels
                for ch = 1:size(Data,1)
                    %Data(ch,ch,:,:) = 0;
                end
                Data = Data./max(Data(:));
            end
            % base-line correction? does it make sense here?
            Databm = repmat(mean(Data(:,:,:,100:376),4),[1 1 1 size(Data,4)]);
            Databs = repmat(std(Data(:,:,:,100:376),[],4),[1 1 1 size(Data,4)]);
            Datap = Data; %(Data-Databm)./Databs;
            plot individual animal resaults
            FIG = dynet_connplot(Datap,tsec,fvec(opt.Freqband(1):opt.Freqband(2)),labels, [], [], [],0);
            set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
            export_fig(FIG,fullfile(opt.figpath,['STOKPDC_' animID{subj} '_normalize_' opt.Normalize]),'-pdf');
            close all;
            
            DataM(:,:,:,:,subj) = Data;%(Data-Databm)./Databs;
        end
        
        % plot average over all animals
        Data = mean(DataM(:,:,:,:,:),5);
        FIG = dynet_connplot(Data./max(Data(:)),linspace(-300,300,size(PDC.mID_40,4)),fvec(opt.Freqband(1):opt.Freqband(2)),labels, [], [], [],0);
        set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
        export_fig(FIG,fullfile(opt.figpath,['STOKPDC_AverageAll_normalize_' opt.Normalize]),'-pdf');
        
        % generate a movie for average FC? (Maybe gamma [50 100])
        
    end
    
end

