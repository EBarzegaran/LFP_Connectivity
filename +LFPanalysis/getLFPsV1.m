function getLFPsV1(LFPpath,infoPath,savePath,Chans)
% retrieve LFP epochs  
% per animal, conditions 
% save 

% INPUT:
    % Chans: ['Layers']/'All';

%% default values

if ~exist('Chans','var') || isempty(Chans) || ~(sum(strcmpi(Chans,{'All','Layers'})))
    Chans = 'Layers';
end

%%

load(infoPath);% load info data
 
% PARAMETERS OF STUDY
animals = unique(cell2mat({dSets.mouse_counter}));
layers = [dSets(1).layer_labels];
conditions = [dSets(1).gratConds.grat_width];

srate=1250;
%fvec = 1:150; % start:end freq, Hz

stf=round(500/(1000/srate)); % start tf before marker
etf=round(500/(1000/srate)); % end tf, after marker
tsec = (-stf : etf) * 1000/srate;

switch Chans
    case 'Layers'
        nchan=numel(layers);
    case 'All'
        nchan=dSets(1).channelCount;
end
       
mc = [dSets.mouse_counter]; % to retrieve nr of exp per mouse

% set 1 Hz high pass filter to the lfp
% [b,a] = butter(2,1/(srate/2),'high');

%x = animals(1); % select animal
for x =animals
    fm=find([dSets.mouse_counter]==x);     % indices the experiments done on this mouse
    animalcnt = find(animals==x);
    disp(['animal ',int2str(animalcnt),': ',int2str(x)])

    r = dSets(fm); % this animal's data, all experiments
    trialPerCond = length(r(1).gratTrials)/numel([r(1).gratConds.grat_width]); % trials per condition recorded

    % init lfp storage for all conditions, exps, one mouse    
    lfpMouse = zeros(stf+etf+1,nchan,trialPerCond,length(r),numel(conditions)); 

    % collect all LFPs, for each condition, each experiment, each recording file
    for exp = 1:numel([r.lfpNames])
        switch Chans
            case 'Layers'
                selElec = cell2mat(r(exp).layer_labels);
            case 'All'
                selElec = 1:r(exp).channelCount;
        end
        fName=char(r(exp).lfpNames);
        disp(['#-------- ' fName ' --------#'])
        lfp = readLFP(fullfile(LFPpath, fName));

        for gNum= 1:numel(conditions) % gratings
            trialNrs=find([r(exp).gratTrials.grat_num]==gNum); % the trials 
            for i = 1:trialPerCond  
                onset=round(r(exp).gratTrials(trialNrs(i)).trial_onset/(1000/srate)*1000); %tf            
                epoch=(lfp(:,(onset-stf):(onset+etf))); % reference timepoint, 32 channel data
                bepoch = CalcBipolar(epoch); % calculate bipolar referenced data
                lfpMouse(:,:,i,exp,gNum) = bepoch(selElec,:)';
            end
        end         
    end
    save(fullfile(savePath,['mID_' num2str(x) '_bipLFPs_' Chans]));
end

%save([outFolder, 'dSets.mat'], 'dSets');
end
