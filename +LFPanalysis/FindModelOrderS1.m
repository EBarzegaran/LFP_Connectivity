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
        'figpath'  ,fullfile(Projfolder,'Results') ...
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

            epochs = permute(LFP.lfpRat(:,:,:), [3,2,1]); %trials, nodes, time
            ff  = 0.99;
            for p =opt.minorder:opt.maxorder
                if mod(p,round((opt.maxorder-opt.minorder)/10))==1, disp([num2str(round((p-opt.minorder)/(opt.maxorder-opt.minorder)*100)) ' %']);end
                REV1(p,1) = LSK_REV1(epochs(1:6,1:6,:), p, ff);    % Y-Y'
                REV2(p,1) = LSK_REV2(epochs(1:6,1:6,:), p, ff);    % using model residuals
                REV1(p,2) = LSK_REV1(epochs(7:12,7:12,:), p, ff);    % Y-Y'
                REV2(p,2) = LSK_REV2(epochs(7:12,7:12,:), p, ff);    % using model residuals
            end
            ModOrds.(animID{subj}).REV1 = REV1;
            ModOrds.(animID{subj}).REV2 = REV2;

        end

        save(fullfile(Projfolder,'ModelOrders.mat'),'ModOrds','opt');
    else
        load(fullfile(Projfolder,'ModelOrders.mat'));
    end
    % add a saving option, since it is slow to estimate the model orders
 %%   Plot the results
 
 ROIs = {'cS1','iS1'};
    if opt.plotfig
        
        for l = 1:2
            Fig = figure;
            for subj = 1:numel(animals)
                subplot(2,3,subj);
                Fields = fields(ModOrds.(animID{subj}));
                % first REV1 and REV2
                sp(1) = plot(opt.minorder:opt.maxorder,Normalize(1-ModOrds.(animID{subj}).REV1(opt.minorder:opt.maxorder,l)),'r');
                hold on; sp(2) = plot(opt.minorder:opt.maxorder,Normalize(1-ModOrds.(animID{subj}).REV2(opt.minorder:opt.maxorder,l)),'b');

               % axis control 
                ylim([-0 1.2]);
                xlim([opt.minorder opt.maxorder]);
                line([opt.minorder opt.maxorder], [0 0],'color','k','linestyle','--')
                set(gca,'yticklabels',[]);

               [~,I1] = max(Normalize(1-ModOrds.(animID{subj}).REV1(:,l)));
               [~,I2] = max(Normalize(1-ModOrds.(animID{subj}).REV2(:,l)));
               xind = opt.minorder:opt.maxorder;
                vline2([xind(I1) xind(I2)],{'--r','--b'});

                title([strrep(animID{subj},'_','-') ]);
                % Plot legends and info on model orders
                if subj==numel(animals)
                    LG = legend(sp,Fields(1:2));
                end  
            end

            set(Fig,'unit','inch','position',[0 0 12 8],'color','w');
            export_fig(Fig,fullfile(opt.figpath,['ModelOrders' '_' ROIs{l}]),'-pdf');
           
        end
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

