function EstimatePDC_STOK(Projfolder,varargin)

%% set default values
    opt = ParseArgs(varargin, ...
        'plotfig'  ,true, ...
        'ModOrds'  ,15,...
        'Cnd'      , 6, ...
        'recalc'   ,false, ...
        'figpath'  ,fullfile(fileparts(fileparts(Projfolder)),'Results') ...
        );

    animals = dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals = {animals.name};
    animID = cellfun(@(x) x(1:6),animals,'uni',false); % animal IDs
    
%% Estimate MVAR parameters using STOK algorithm
    xticklabels = -500:100:500;
    xticks = linspace (0,1250,numel(xticklabels));
    
    if ~exist(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd)]),'file') || opt.recalc
        for subj = 1:numel(animals)
            
            %---------------------- Load Data -----------------------------
            disp(animID{subj});
            LFP = load(fullfile(Projfolder,animals{subj}),'lfpMouse');
            labels = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
            dims = size(LFP.lfpMouse);
            LFP.lfpMouse = reshape(LFP.lfpMouse, dims(1),dims(2),dims(3)*dims(4),dims(5)); %treat all exps (recordings) as one

            %--------------- Estimate MVAR params and PDC -----------------
            Cond = opt.Cnd;
            % Prepare the data
            epochs = permute(LFP.lfpMouse(:,:,:,Cond), [3,2,1]); %trials, nodes, time
            % STOC and PDC
            ff=.98;
            KF = dynet_SSM_SALKff(epochs(:,:,250:1001),opt.ModOrds,ff);  
            PDC = dynet_ar2pdc(KF,srate,fvec,measure,keepdiag,flow);
        end

        save(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd)]),'');
    else
        load(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd)]));
    end

end