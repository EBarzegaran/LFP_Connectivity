function [ModOrds] = FindModelOrderS1(Projfolder,varargin)
% determine optimal p for Kalman based state space modeling using saved data, no high-pass
% Author: Plomp
    % Latest modification: Elham Barzegaran

%% set default values
    opt = ParseArgs(varargin, ...
        'plotfig'  ,true, ...
        'maxorder' ,30,...
        'minorder' ,1,...
        'recalc'     ,true, ...
        'figpath'  ,fullfile(Projfolder,'Results'),...
        'TimeWin'    ,[-100 100] ...
        );

    animals = dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals = {animals.name};
    animID = cellfun(@(x) x(1:4),animals,'uni',false); % animal IDs
    
%% Compute model order for each animal
    xticklabels = -100:50:100;
    xticks = linspace (0,1250,numel(xticklabels));
    disp('Estimating optimal model orders:');
    
    if ~exist(fullfile(Projfolder,'ModelOrders.mat'),'file') || opt.recalc
        for subj = 1:numel(animals)
            disp(animID{subj});
            LFP = load(fullfile(Projfolder,animals{subj}),'lfpRat');
            labels = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
            dims = size(LFP.lfpRat);
            
            %---------------------------- optimize model order-----------------
            % QUESTION: which part of data should be used? i.e. which time window

            epochs	= permute(LFP.lfpRat(:,:,:), [3,2,1]); %trials, nodes, time
            freqs	=   1:150;
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
            ff = .99;
            if size(freqs,2)>size(freqs,1); freqs=freqs'; end
            for p = opt.minorder:opt.maxorder
                disp(['Model Order: ' num2str(p)]);
                KF       = dynet_SSM_STOK(epochs(:,:,tvec),p,ff);
                Z        = repmat(exp(-2*pi*1i*freqs/srate),[1 p]).^repmat((1:p),[numel(freqs) 1]);
                A        = repmat(eye(nodes), [1 1 length(freqs) Time]);
                
                for k = 1:p
                    tmp  = repmat(-KF.AR(:, :, k,:), [1 1 length(freqs) 1]);
                    A    = A + bsxfun(@times,tmp,reshape(Z(:,k),1,1,[]));
                end
                
                NCov = 1;
                
                PSD_p{p} = (1/srate*NCov)./abs(repmat(eye(size(A,1)),1,1,size(A,3),size(A,4))+A).^2;
                
            end
            PSD_param{subj} = PSD_p{40}(:,:,20:end,:);
            
            
            %-------------------------Estimate the MSE of PSDs------------------
            SN = abs(S(:,:,5:end,:))./(sum(sum(abs(S(:,:,5:end,:)),3),4));
            
            PSD_nonparam{subj} = abs(S(:,:,20:end,:));
            
            %FIG = figure;
            for p = opt.minorder:opt.maxorder
                PSD         = PSD_p{p}(:,:,5:end,:);
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
        save(fullfile(Projfolder,['ModelOrders' '_PSD.mat']),'MSE','freqs','opt','tsec');
        tsec1 = tsec(200:1000);
        freqs1 = freqs(20:end);
        save(fullfile(Projfolder,'PSD_param_nonpram.mat'),'PSD_nonparam','PSD_param','tsec1','freqs1');
    else
        load(fullfile(Projfolder,['ModelOrders' '_PSD.mat']));
        load(fullfile(Projfolder,'PSD_param_nonpram.mat'));
    end
    % add a saving option, since it is slow to estimate the model orders

%% Plot the results
    FIG2 = figure;
    for subj = 1:numel(animals)
        subplot(2,3,subj),
        plot(opt.minorder:opt.maxorder,MSE(opt.minorder:opt.maxorder,subj));
        [~,I] = min(MSE(:,subj));
        vline2(I,{'--k'});
        title([strrep(animID{subj},'_','-') ' - Order = ' num2str(I)]);
    end 
    set(FIG2,'unit','inch','position',[0 0 20 10],'color','w');
    export_fig(FIG2,fullfile(opt.figpath,['ModelOrders_' '_PSD']),'-pdf'); 

    %%
    PSD_paramM = mean(cat(5,PSD_param{1:6}),5);
    PSD_nonparamM = mean(cat(5,PSD_nonparam{1:6}),5);
    labels = arrayfun(@(x) ['L' num2str(x)],1:12,'uni',false);
    %%
    H = {'cS1','iS1'};
    for l = 1:2
    FIG3 = figure;
        for c = (6*(l-1)+1):6*l
            
            subplot(6,2,(c-((l-1)*6))*2-1)
            imagesc(abs(squeeze(PSD_nonparamM(c,c,:,:))));
            set(gca,'xtick',2:100:numel(tsec),'xticklabel',round(tsec(2:100:end)))
            vline(626,'w--')
            axis xy;
            title(labels{c})

            subplot(6,2,(c-((l-1)*6))*2)
            imagesc(abs(squeeze(PSD_paramM(c,c,:,:))));
            set(gca,'xtick',1:100:numel(tsec),'xticklabel',round(tsec(2:100:end)))
            vline(626,'w--')
            axis xy;
            title(labels{c})

        end
        set(FIG3,'unit','inch','position',[2 2 10 16],'color','w')
        export_fig(FIG3,fullfile(opt.figpath,['PSD_ModelOrders10_' '_PSD_' H{l}]),'-pdf');
    end
end

