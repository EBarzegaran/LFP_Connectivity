
function visualizeLFP(Projfolder,varargin)
% This function visualize mouse V1 LFP data, individuals and average, with
% no prior preprocessing, i.e. filtering
% Syntax: visualizeLFP(Projfolder,varargin)

% INPUTs:
    % Projfolder: A string containing a string with folder path where preprocessed lfp files are stored

%<OPTIONS>:
    % 'savefig': Indicates if the figures should be saved or not, [true]/false
    % 'figpath': A string containing a string with folder path where the figures will be saved
    
% Author: Elham Barzegaran, 16/07/2019    
   
    %% Set default values
    opt = ParseArgs(varargin, ...
        'savefig'  ,'true', ...
        'figpath'  ,fullfile(fileparts(fileparts(Projfolder)),'Results') ...
        );

    animals = dir(fullfile(Projfolder,'*Layers.mat'));
    animals = {animals.name};
    animID = cellfun(@(x) x(1:6),animals,'uni',false);
    
    
    %% Plot individual and average LFPs
    
    Colors = brewermap(6,'Spectral');
    offsets =.03*5:-.03:0;
    xticklabels = -500:100:500;
    xticks = linspace (0,1250,numel(xticklabels));
    Fig = figure;

    % Plot individual animal data
    for subj = 1:numel(animals)
        LFP.(animID{subj}) = load(fullfile(Projfolder,animals{subj}),'lfpMouse');
        subplot(3,3,subj)
        for l = 1:6
            hold on;plot(squeeze(mean(mean(LFP.(animID{subj}).lfpMouse(:,l,:,:,6),3),4))+offsets(l),'Color',Colors(l,:))
        end
        Mdata(:,:,subj) = squeeze(mean(mean(LFP.(animID{subj}).lfpMouse(:,:,:,:,6),3),4));
        nTrl(subj) = size(LFP.(animID{subj}).lfpMouse,3) * size(LFP.(animID{subj}).lfpMouse,4);
        title([strrep(animID{subj},'_','-') ', nTr = ' num2str(nTrl(subj))])
        xlim([250 1001]);
        set(gca,'xtick',xticks,'xticklabel',xticklabels,'yticklabel',[]);
        line([626 626],[-.05 .5],'linestyle','-.','linewidth',1.2,'color','k')
        ylim([-.05 .21])

    end

    % Plot average LFP over 7 animals
    subplot(3,3,8)
    for l = 1:6
       hold on;plot(squeeze(mean(Mdata(:,l,:),3))+offsets(l),'Color',Colors(l,:))
    end
    xlim([250 1001]);
    set(gca,'xtick',xticks,'xticklabel',xticklabels,'yticklabel',[]);
    line([626 626],[-.05 .5],'linestyle','-.','linewidth',1.2,'color','k')
    ylim([-.05 .2])
    title('Average'); xlabel('Time(s)', 'fontweight','bold')
    set(gcf,'unit','inch','position',[0 0 20 15],'color','w')
    legend(arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));

    if opt.savefig
        if ~exist(opt.figpath,'dir')
            try
                mkdir(opt.figpath)
            catch
                error('Invalid path for saving figures')
            end
        end
        export_fig(Fig,fullfile(opt.figpath,'LFP-nonfiltered'),'-pdf');
    end

end

