function getLFPsS1(LFPpath,savePath,Chans)
% retrieve LFP epochs from S1 Rat dataset
% per animal, conditions 
% save 

% INPUT:
    % Chans: 'All';

%% default values

if ~exist('Chans','var') || isempty(Chans) || ~(sum(strcmpi(Chans,{'All','Layers'})))
    Chans = 'All';
end

if ~exist(savePath,'dir')
    mkdir(savePath)
end
%%
filenames = dir(fullfile(LFPpath,'*.SEF'));
animals = {filenames.name};
for x = 1:numel(animals)
    sefFile=fullfile(LFPpath,animals{x});
    [nchan,nbsamples,EEG.fSamp,orgDat] = opensef(sefFile);
    nchan = 16;
    orgDat=orgDat(:,[1:4,6:12,14:16]);% remove trigger channel, plus one channel on the right to get equal nrs
    nchan=nchan-2;
    orgDat = orgDat(:,[2:7, 9:14]); % exclude pia electrode, leaving bilateral L1-L6
    nchan=nchan-2;

    % load condition markers
    % N.B. mrk tfs are zero based, whereas matlab is 1-based
    mrkFile=[sefFile '.mrk'];
    fid = fopen(mrkFile);
    header= textscan(fid,'%s', 1);
    mrks=textscan(fid,'%d64 %d64 %s');
    fclose(fid)


    srate = EEG.fSamp;
    pdc_stf=100/(1000/srate); % start tf before marker, 100 ms 
    pdc_etf=100/(1000/srate); % end tf, after marker, 100 ms

    %----------------------------------- select data
    
    nTrials=size(mrks{1},1);
    epochs=zeros(nTrials, nchan, pdc_stf+pdc_etf+1); % Y = data ( trials, channels, samples)
    for d = 1:nTrials % nr of trials, left/right
        mrkTf=mrks{1}(d)+1; % tf of this marker
        dat=orgDat([(mrkTf-pdc_stf):(mrkTf+pdc_etf)],:)';

        if strcmp(mrks{3}(d),'"Left"') %swap left/right so 1:6 is contralateral S1 layers
            dat = dat([7:12, 1:6], :);
        end

        epochs(d,:,:) = dat;        
        disp(['Collect epochs , did  ', int2str(d) ' out of ' num2str(nTrials)])
    end
    %     % plot epochs
    %     avEp=squeeze(mean(epochs,1));
    %     tsec = (-pdc_stf:(size(epochs,3)-pdc_stf-1)) * (1000/srate);
    %     plot(tsec,avEp')
    %     legend('location', 'southeast')
    %     
    % ------------------------- bipolar rerefrencing ----------------------
    
    % WHY NO BIPOLAR REFERENCING HERE?
%      bepoch(1:6,:,:) = CalcCSD(permute(epochs(:,1:6,:),[2 3 1])); % calculate bipolar referenced data
%      bepoch(7:12,:,:) = CalcCSD(permute(epochs(:,7:12,:),[2 3 1])); % calculate bipolar referenced data
%     bepochs(:,1:6,:)  = epochs(:,1:6,:) - mean(epochs(:,1:6,:),2);
%     bepochs(:,7:12,:)  = epochs(:,7:12,:) - mean(epochs(:,7:12,:),2);
%     lfpRat = permute(epochs,[3 2 1]);
    tsec = (-pdc_stf:(size(epochs,3)-pdc_stf-1)) * (1000/srate);

    save(fullfile(savePath,[animals{x}(1:4) '_LFPs_' Chans]),'lfpRat','srate','nTrials','nchan','tsec');
end

end
