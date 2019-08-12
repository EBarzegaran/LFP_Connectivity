% Main script for analysis of mouse V1 LFP data, for more info on the
% dataset please see: Plomp, et al. 2019, JNeuro

% Author: Elham Barzegaran, July, 2019
        % e.barzegaran@gmail.com
%% Add toolboxes
clear; clc;
mac = false;
if mac
    addpath(genpath('/Users/elhamb/switchdrive/Institution/Codes'));
    addpath(genpath('/Users/elhamb/Documents/Codes/Git/EEGSourceSim'));
    MainDir = '/Users/elhamb/switchdrive/Institution/Data/mouseV1_prep';
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
LFPanalysis.visualizeLFP(Projfolder,    'savefig',true, 'figpath',fullfile(Projfolder,'Results'));

%% Run parametric and non-parametric connectivity analysis and compare the results

% (1) Find optimal model order for parametric methods: I start with stoc algorithm
ModOrds = LFPanalysis.FindModelOrderV1(Projfolder,'recalc',true);
%% Estimate PDC, calculate statistics based on permutation test
% (2) estimate pdc based on MVAR coefficients
LFPanalysis.EstimatePDC_STOK(Projfolder,    'Normalize','Channel',  'recalc',true,  'plotfig',false);

%% other analysis: all channels-> decomposition like SSD, RCA and etc.
% decide later if you want to implement this
