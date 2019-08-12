
function visualizeLFPS1(Projfolder,varargin)
% This function visualize Rat S1 LFP data, individuals and average, with
% no prior preprocessing, i.e. filtering
% Syntax: visualizeLFP(Projfolder,varargin)

% INPUTs:
    % Projfolder: A string containing a string with folder path where preprocessed lfp files are stored

%<OPTIONS>:
    % 'savefig': Indicates if the figures should be saved or not, [true]/false
    % 'figpath': A string containing a string with folder path where the figures will be saved
    
% Author: Elham Barzegaran, 12/08/2019    
   
    %% Set default values
    opt = ParseArgs(varargin, ...
        'savefig'  ,'true', ...
        'figpath'  ,fullfile(fileparts(fileparts(Projfolder)),'Results') ...
        );

    animals = dir(fullfile(Projfolder,'*Layers.mat'));
    animals = {animals.name};
    animID = cellfun(@(x) x(1:4),animals,'uni',false);
    
    
    %% Plot individual and average LFPs
    
    Colors = brewermap(6,'Spectral');
    offsets =[0:1000:5000 7000:1000:12000];
    xticklabels = -100:50:100;
    xticks = linspace (0,1250,numel(xticklabels));
    Fig = figure;
    Labels = [arrayfun(@(x) ['cS1-L' num2str(x)],1:6,'uni',false),arrayfun(@(x) ['iS1-L' num2str(x)],1:6,'uni',false)];
    % Plot individual animal data
    for subj = 1:numel(animals)
        LFP.(animID{subj}) = load(fullfile(Projfolder,animals{subj}),'lfpRat');
        subplot(2,3,subj)
        for l = 1:6
            hold on;plot(squeeze(mean(LFP.(animID{subj}).lfpRat(:,l,:),3))+offsets(l),'Color',Colors(l,:))
        end
        
        for l = 1:6
            hold on;plot(squeeze(mean(LFP.(animID{subj}).lfpRat(:,l+6,:),3))+offsets(l+6),'Color',Colors(l,:))
        end
        %Mdata(:,:,subj) = squeeze(mean(mean(LFP.(animID{subj}).lfpMouse(:,:,:,:,6),3),4));
        nTrl(subj) = size(LFP.(animID{subj}).lfpRat,3) * size(LFP.(animID{subj}).lfpRat,4);
        title([strrep(animID{subj},'_','-') ', nTr = ' num2str(nTrl(subj))])
        xlim([0 1251]);
        set(gca,'xtick',xticks,'xticklabel',xticklabels,'ytick',offsets,'yticklabels',Labels);
        line([626 626],[-3000 15000],'linestyle','-.','linewidth',1.2,'color','k')
        ylim([-3000 15000])

    end

    % Plot average LFP over 7 animals
     xlabel('Time(s)', 'fontweight','bold')
    set(gcf,'unit','inch','position',[0 0 20 15],'color','w')
    %legend();

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

