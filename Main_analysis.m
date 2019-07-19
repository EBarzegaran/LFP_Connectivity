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
else
    addpath(genpath('E:\Elham\Codes\Git\EEGSourceSim')); % I use some functions from here
    addpath(genpath('E:\Elham\Codes\NonGit\mvgc_v1.0'));
    addpath(genpath('C:\Users\Barzegar\switchdrive\Institution\Codes\dynet_toolbox-master-'));
end

% Set the paths
MainDir = fileparts(mfilename('fullpath'));
LFPpath = fullfile(MainDir,'data','LFPs'); % path of the original LFP files
Projfolder = fullfile(MainDir,'data','bipLFPs');% path to LFP bipolar data
infofile = fullfile(MainDir,'data','dInfo','dInfoCombined.mat'); % session info files
addpath(genpath(MainDir));

%% extract the full dataset (all channels included)
LFPanalysis.getLFPs(LFPpath,infofile,Projfolder,'Layers');

%% visualize LFP data
LFPanalysis.visualizeLFP(Projfolder,'savefig',true,'figpath',fullfile(Projfolder,'Results'));

%% other analysis: all channels-> decomposition like SSD, RCA and etc.
% decide later if you want to implement this


%% Run parametric and non-parametric connectivity analysis and compare the results

% (1) Find optimal model order for parametric methods: I start with stoc algorithm
ModOrds = LFPanalysis.FindModelOrder(Projfolder,'recalc',true);
%%
% (2) estimate pdc based on MVAR coefficients
LFPanalysis.EstimatePDC_STOK(Projfolder);
%%

