function [ModOrds] = FindModelOrderV1(Projfolder,varargin)
% determine optimal p for Kalman based state space modeling using saved data, no high-pass
% Author: Plomp
    % Latest modification: Elham Barzegaran

%% set default values
    opt = ParseArgs(varargin, ...
        'plotfig'  ,true, ...
        'maxorder' ,30,...
        'Cnd'      , 6, ...
        'recalc'     ,false, ...
        'figpath'  ,fullfile(fileparts(fileparts(Projfolder)),'Results') ...
        );

    %animals = dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals = dir(fullfile(Projfolder,'mID*Layers.mat')); % files in the project folder
    animals = {animals.name};
    animID = cellfun(@(x) x(1:6),animals,'uni',false); % animal IDs
    
    if ~exist(opt.figpath)
        mkdir(opt.figpath);
    end
    
%% Compute model order for each animal
    disp('Estimating optimal model orders:');
    
    if ~exist(fullfile(Projfolder,['ModelOrders_Cnd' num2str(opt.Cnd) '.mat'] ),'file') || opt.recalc
        for subj = 1:numel(animals)
            disp(animID{subj});
            LFP = load(fullfile(Projfolder,animals{subj}),'lfpMouse');
            labels = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
            dims = size(LFP.lfpMouse);
            LFP.lfpMouse = reshape(LFP.lfpMouse, dims(1),dims(2),dims(3)*dims(4),dims(5)); %treat all exps (recordings) as one

            %---------------------------- optimize model order-----------------
            % QUESTION: which part of data should be used? i.e. which time window
            Cond = opt.Cnd;
            epochs = permute(LFP.lfpMouse(:,:,:,Cond), [3,2,1]); %trials, nodes, time
            %save(fullfile(Projfolder,'temp1.mat'),'epochs')
%             [aic,bic,moaic,mobic] = tsdata_to_infocrit(permute(epochs(:,:,250:1001),[2 1 3]),opt.maxorder,'LWR',false); % stationary
            [aic,bic,moaic,mobic] = tsdata_to_infocrit(permute(epochs(:,:,1:625),[2 1 3]),opt.maxorder,'LWR',false); % stationary, on prestimulus data
            ff  = 0.99;
            for p = 1:opt.maxorder
                if mod(p,3)==0, disp([num2str(round(p/opt.maxorder*100)) ' %']);end
                REV1(p) = LSK_REV1(epochs(:,:,:), p, ff);    % Y-Y'
                REV2(p) = LSK_REV2(epochs(:,:,:), p, ff);    % using model residuals
            end
            ModOrds.(animID{subj}).REV1 = REV1;
            ModOrds.(animID{subj}).REV2 = REV2;
            ModOrds.(animID{subj}).aic = aic;
            ModOrds.(animID{subj}).bic = bic;
            ModOrds.(animID{subj}).Cnd = opt.Cnd;
        end

        save(fullfile(Projfolder,['ModelOrders_Cnd' num2str(opt.Cnd)]),'ModOrds');
    else
        load(fullfile(Projfolder,['ModelOrders_Cnd' num2str(opt.Cnd)]));
    end
    % add a saving option, since it is slow to estimate the model orders
 %%   Plot the results
    if opt.plotfig
        Fig = figure;
        for subj = 1:numel(animals)
            subplot(3,3,subj);
            Fields = fields(ModOrds.(animID{subj}));
            % first REV1 and REV2
            sp(1) = plot(Normalize(1-ModOrds.(animID{subj}).REV1),'r');
            hold on; sp(2) = plot(Normalize(1-ModOrds.(animID{subj}).REV2),'b');
            % plot aic and bic
            sp(3) = plot(Normalize(ModOrds.(animID{subj}).aic),'g');
           sp(4) = plot(Normalize(ModOrds.(animID{subj}).bic),'k');
           
           % axis control 
            ylim([-1.2 1.2]);
            xlim([1 opt.maxorder]);
            line([1 opt.maxorder], [0 0],'color','k','linestyle','--')
            set(gca,'yticklabels',[]);
            
           [~,I1] = min(Normalize(ModOrds.(animID{subj}).aic));
           [~,I2] = min(Normalize(ModOrds.(animID{subj}).bic));
            vline2([I1 I2],{'--g','--k'});
            
           [~,I1] = max((1-ModOrds.(animID{subj}).REV1));
           [~,I2] = max((1-ModOrds.(animID{subj}).REV2));
            vline2([I1 I2],{'--r','--b'});
            
            title([strrep(animID{subj},'_','-') ' - Cond #' num2str(opt.Cnd)]);
            % Plot legends and info on model orders
            if subj==numel(animals)
                LG = legend(sp,Fields(1:4));
            end  
        end
        
        set(Fig,'unit','inch','position',[0 0 20 15],'color','w');
        export_fig(Fig,fullfile(opt.figpath,['ModelOrders_' num2str(opt.Cnd)]),'-pdf');
    end
    

end

function normdata = Normalize(data)
% Normalizing a vector purly for visualization purpose! Do not use this
% function in computations.
% data should be a vectror not a matrix

if sum(data>=0)==numel(data)
    normdata = (data-min(data))/(max(data)-min(data));
end

if sum(data<=0)==numel(data)
    pdata = abs(data);
    normdata = -(pdata-min(pdata))/(max(pdata)-min(pdata));
end
end

