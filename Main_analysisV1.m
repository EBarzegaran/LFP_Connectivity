% Main script for analysis of mouse V1 LFP data, for more info on the
% dataset please see: Plomp, et al. 2019, JNeuro

% Author: Elham Barzegaran, July, 2019
        % e.barzegaran@gmail.com
%% Add toolboxes
clear; clc;
mac = true;
if mac
    addpath(genpath('/Users/elhamb/switchdrive/Institution/Codes'));
    addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
    MainDir = '/Users/elhamb/switchdrive/Institution/Data/mouseV1_prep';
    addpath(genpath('/Users/elhamb/Documents/Codes/NonGit/nonparametricGGC_toolbox'));
else
    addpath(genpath('E:\Elham\Codes\Git\EEGSourceSim')); % I use some functions from here
    addpath(genpath('E:\Elham\Codes\NonGit\mvgc_v1.0'));
    addpath(genpath('C:\Users\Barzegar\switchdrive\Institution\Codes\dynet_toolbox-master-'));
    MainDir = 'E:\Experiments\mouseV1\';
end

% Set the paths
LFPpath = fullfile(MainDir,'data','LFPs'); % path of the original LFP files
Projfolder = fullfile(MainDir,'data','bipLFPs');% path to LFP bipolar data
infofile = fullfile(MainDir,'data','dInfo','dInfoCombined.mat'); % session info files
addpath(genpath(fileparts(mfilename('fullpath'))));

%% extract the full dataset (all channels included)
LFPanalysis.getLFPsV1(LFPpath,infofile,Projfolder,'Layers');

%% visualize LFP data
LFPanalysis.visualizeLFP(Projfolder,...
    'savefig',true,...
    'figpath',fullfile(Projfolder,'Results'));

%% Run parametric and non-parametric connectivity analysis and compare the results

% (1) Find optimal model order for parametric methods: I start with stok algorithm
ModOrds = LFPanalysis.FindModelOrderV1(Projfolder,...%Projfolder,...
    'recalc',true);

%% (1.5) Find the optimal model order by comparing parametric and
% non-paramnetric PSDs
ModOrds2 = LFPanalysis.FindModelOrderV1PSD(Projfolder,...
    'recalc'    ,true,...
    'maxorder'  ,20);

%% Estimate PDC, calculate statistics based on permutation test
% (2) estimate pdc based on MVAR coefficients
for cnd = 1:10
    disp(['Condition: ' num2str(cnd)])
    [~, Diff{cnd},fvec,tsec] = LFPanalysis.EstimatePDC_STOKV1(Projfolder,...
        'Normalize'     ,'Channel',...
        'Cnd'           ,cnd,...
        'NormPrestim'   ,true,...
        'recalc'        ,true,...
        'plotfig'       ,false,...
        'ModOrds'       ,15,...
        'IDs'           ,[],...
        'TimeWin'       ,[-200 300],...
        'Freqband'      ,[20 150],...
        'doPermute'     ,false);
    close all;
end


%%

Diffmat = cellfun(@(x) cat(3,x{:}),Diff,'uni',false);
FIG = figure;
CS = brewermap(10,'YlOrRd');
for l = 2:5
    subplot(4,1,l-1);
    hold on;
    line([0 626],[0 0],'linestyle','--','color','k')
    
    for cnd = 2:2:numel(Diffmat)
        Data = (squeeze(mean(Diffmat{cnd}(:,:,l),1)));
        lg(cnd/2)= plot(tsec,Data,'color',[CS(cnd,:) 1],'linewidth',1.5);
        
    end
    if l==5
        legend(lg,arrayfun(@(x) ['Cnd-' num2str(x)],2:2:10,'uni',false))
          xlabel('Time(msec)');
        ylabel('Directionality (Upward+)')
    end
    xlim([-100 300])
    ylim([-.1 0.1])
    title(['L' num2str(l)])
    vline(0,'k--');
    set(gca,'fontsize',12)
end

set(FIG,'unit','inch','position',[1 1 15 15],'color','w');
export_fig(FIG,fullfile(MainDir,'Results','Directionality_PDC_allcond'),'-pdf')
%%
ind = (tsec>=-100);
tseccorr = tsec(ind);
xtickso = find(ismember(tseccorr,-100:100:300));
ytickso = find(ismember(fvec,20:40:150));


locs = [1:2:10 2:2:10];
layers = {2,3,4,5,2:5};
for l =5:numel(layers)
    FIG2 = figure;
    for cnd = 1:10
        subplot(5,2,locs(cnd));
        imagesc(mean(Diffmat{cnd}(:,ind,layers{l}),3));
        axis xy
        caxis([-.07 .07])
        set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',ytickso,'yticklabel',20:40:150);
        vline(find(tseccorr==0),'k--')    
        if cnd==5
           xlabel('time(msec)');
           ylabel('Freq(Hz)')
        end
        title(['Cnd-' num2str(cnd)])
    end
    colormap(jmaColors('coolhotcortex'));
    set(FIG2,'unit','inch','position',[1 1 15 15],'color','w');
end

export_fig(FIG2,fullfile(MainDir,'Results','Directionality_PDC_allcond_allfreqs'),'-pdf')
%% other analysis: all channels-> Decompositions like SSD, RCA and etc.
% decide later if you want to implement this

