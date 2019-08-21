function ModOrds = FindModelOrderV1PSD(Projfolder,varargin)

%% Reference for non-parametric PSD estimations:
% [1] Baccalá, et al. Partial directed coherence: a new concept in neural
    % structure determination, 2001


%% set default values
    opt = ParseArgs(varargin, ...
        'plotfig'   ,true, ...
        'maxorder'  ,30,...
        'Cnd'       , 6, ...
        'recalc'    ,false, ...
        'TimeWin'   ,[],...
        'figpath'   ,fullfile(fileparts(fileparts(Projfolder)),'Results') ...
        );

    %animals = dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals = dir(fullfile(Projfolder,'mID*Layers.mat')); % files in the project folder
    animals = {animals.name};
    animID  = cellfun(@(x) x(1:6),animals,'uni',false); % animal IDs
    
    if ~exist(opt.figpath)
        mkdir(opt.figpath);
    end
    
    ModOrds = [];
%% Compute model order for each animal
    disp('Estimating optimal model order:');
    
    if ~exist(fullfile(Projfolder,['ModelOrders_Cnd' num2str(opt.Cnd) '_PSD.mat'] ),'file') || opt.recalc
        for subj = 1:numel(animals)
            disp(['Animal ID: ' animID{subj}]);
            LFP          =  load(fullfile(Projfolder,animals{subj}),'lfpMouse');
            labels       =  arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
            dims         =  size(LFP.lfpMouse);
            LFP.lfpMouse =  reshape(LFP.lfpMouse, dims(1),dims(2),dims(3)*dims(4),dims(5)); %treat all exps (recordings) as one

            %---------------------------- optimize model order-----------------
            % QUESTION: which part of data should be used? i.e. which time window
            Cond                    =   opt.Cnd;
            epochs                  =   permute(LFP.lfpMouse(:,:,:,Cond), [3,2,1]); %trials, channels, time
            %[aic,bic,moaic,mobic]   =   tsdata_to_infocrit(permute(epochs(:,:,1:625),[2 1 3]),opt.maxorder,'LWR',false); % stationary, on prestimulus data
            ff                      =   0.99;
            freqs                    =   1:150;
            load(fullfile(Projfolder,animals{subj}),'srate');
            load(fullfile(Projfolder,animals{subj}),'tsec');
            if isempty(opt.TimeWin)
                opt.TimeWin = [min(tsec) max(tsec)];
            end
            tvec                    =   (tsec>=opt.TimeWin(1)) & (tsec<=opt.TimeWin(2));
            tsec                    =   tsec(tvec);
            
            %-----------------(1) estimate non-parametric PSD using wavelet----------
            % Time x Trial x Channel
            [S, ~, ~]  = xwt_cmorl_nv(permute(epochs,[3 1 2]), srate, freqs, 1, 3);
            
            %------------(2) estimate parametric PSD using STOK algorithm -------------
            nodes   = size(epochs,2);
            Time    = size(epochs,3);
            if size(freqs,2)>size(freqs,1); freqs=freqs'; end
            for p = 1:opt.maxorder
                disp(['Model Order: ' num2str(p)]);
                KF       = dynet_SSM_STOK(epochs(:,:,tvec),p,ff);
                Z        = repmat(exp(-2*pi*1i*freqs/srate),[1 p]).^repmat((1:p),[numel(freqs) 1]);
                A        = repmat(eye(nodes), [1 1 length(freqs) Time]);
                
                for k = 1:p
                    tmp  = repmat(-KF.AR(:, :, k,:), [1 1 length(freqs) 1]);
                    A    = A + bsxfun(@times,tmp,reshape(Z(:,k),1,1,[]));
                end
                
                % A = bsxfun(@rdivide,(A), sqrt(sum(abs(A).^2,2))); % normalize?
                
%                 PYe = ((KF.PY-epochs(:,:,tvec)).^2); % error covariance?
%                 MSE = mean(PYe(:));
                NCov = 1;
                
                PSD_p{p} = (1/srate*NCov)./abs(repmat(eye(size(A,1)),1,1,size(A,3),size(A,4))+A).^2;
                
            end
            PSD_param{subj} = PSD_p{15}(:,:,20:end,200:1000);
            
            
            %-------------------------Estimate the MSE of PSDs------------------
            SN = abs(S(:,:,5:end,200:1000))./(sum(sum(abs(S(:,:,5:end,200:1000)),3),4));
            
            PSD_nonparam{subj} = abs(S(:,:,20:end,200:1000));
            
            %FIG = figure;
            for p = 1:opt.maxorder
                PSD         = PSD_p{p}(:,:,5:end,200:1000);
                PSDN        = PSD./(sum(sum(PSD,3),4));
                SE          = (PSDN-SN).^2;
                for c = 1:size(SE,1)
                    mse(c) = mean(mean(SE(c,c,:,:)));
                end
                MSE(p,subj) = mean(mse);
%                 subplot(5,6,p),
%                 imagesc(squeeze(PSD_p{p}(3,2,5:end,200:1000)))
%                 title(['p=' num2str(p) ' MSE = ' num2str(round(MSE(p,subj),2))]);
            end
            
        end
        save(fullfile(Projfolder,['ModelOrders_Cnd' num2str(opt.Cnd) '_PSD.mat']),'MSE','freqs','opt','tsec');
        tsec1 = tsec(200:1000);
        freqs1 = freqs(20:end);
        save(fullfile(Projfolder,'PSD_param_nonpram.mat'),'PSD_nonparam','PSD_param','tsec1','freqs1');
    else
        load(fullfile(Projfolder,['ModelOrders_Cnd' num2str(opt.Cnd) '_PSD.mat']));
    end
    %% Plot the results
    FIG2 = figure;
    for subj = 1:numel(animals)
        subplot(3,3,subj),
        plot(MSE(:,subj));
        [~,I] = min(MSE(:,subj));
        vline2(I,{'--k'});
        title([strrep(animID{subj},'_','-') ' - Order = ' num2str(I)]);
    end 
    set(FIG2,'unit','inch','position',[0 0 20 15],'color','w');
    export_fig(FIG2,fullfile(opt.figpath,['ModelOrders_' num2str(opt.Cnd) '_PSD']),'-pdf');
end