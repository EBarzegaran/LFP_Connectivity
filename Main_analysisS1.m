% Main script for analysis of Rat S1 LFP data, for more info on the
% dataset please see: Plomp, et al. 2014, EJN

% Author: Elham Barzegaran, August, 2019
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
    addpath(genpath('C:\Users\Barzegar\switchdrive\Institution\Codes\'));
    MainDir = 'E:\Experiments\RatS1';
end

% Set the paths
LFPpath = fullfile(MainDir,'dataSets'); % path of the original LFP files
Projfolder = fullfile(MainDir,'LFPs');% path to LFP bipolar data
% infofile = fullfile(MainDir,'data','dInfo','dInfoCombined.mat'); % session info files
addpath(genpath(fileparts(mfilename('fullpath'))));

%% extract the full dataset (all channels included)
LFPanalysis.getLFPsS1(LFPpath,Projfolder,'Layers');

%% visualize LFP data
LFPanalysis.visualizeLFPS1(Projfolder,...
    'savefig',true,...
    'figpath',fullfile(Projfolder,'Results'));

%% other analysis: all channels-> decomposition like SSD, RCA and etc.
% decide later if you want to implement this


%% Run parametric and non-parametric connectivity analysis and compare the results
% (1) Find optimal model order for parametric methods: I start with stoc algorithm
ModOrds = LFPanalysis.FindModelOrderS1(Projfolder,...
    'recalc',false,...
    'maxorder',120,...
    'minorder',1);

%% Estimate PDC, calculate statistics based on permutation test
% (2) estimate pdc based on MVAR coefficients
LFPanalysis.EstimatePDC_STOKS1(Projfolder,...
    'Normalize','All',...
    'recalc',false,...
    'plotfig',false,...
    'SStat','Paired',...
    'ModOrds',50,...
    'dopermute',false);

%%

