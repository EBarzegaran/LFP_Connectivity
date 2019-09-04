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
        'recalc'        ,false,...
        'plotfig'       ,false,...
        'ModOrds'       ,15,...
        'IDs'           ,[],...
        'TimeWin'       ,[-200 300],...
        'Freqband'      ,[20 150],...
        'doPermute'     ,false);
    close all;
end

LFPanalysis.plotResultsDirect(Diff,tsec,fvec,MainDir);
%%


%% other analysis: all channels-> Decompositions like SSD, RCA and etc.
% decide later if you want to implement this

