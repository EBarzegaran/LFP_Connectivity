function  [PDC_Results, Direction_Stats]  =   bootstrap_PDC  (epochs,  srate,  fvec,   tsec, labels, ROIs,  varargin)

% This function calculates PDC based on STOK algorithm, and using
% bootstrapping it estimnates the significant results.
% ********** Later I might add cluster-based permutation test to this analysis
%
% Syntax: [PDC_Results, Direction_Stats]  =   bootstrap_PDC  (epochs,  srate,  fvec,   tsec, labels, ROIs,  varargin)
% -------------------------------------------------------------------------
% INPUTS
% - epochs:     A cell array (size = number of subjects) of matrices, containing trial x channel x  time matrix
% - srate:      Sampling rate
% - fvec:       A vector containing the frequency indices of interest
% - tsec:       A vector containing the time points - stimulus locked
%
% - <OPTIONS>
% - 'nboot':        number of bootstraps
% - 'bootsize':     size of bootstrap draws
% - 'ModOrds':      model order of STOK
% - 'ff':           Variance to keep in STOK
% - 'measure':      options: 'PDCnn' , 'PDC', 'sPDC', for more info see function: dynet_ar2pdc 
% - 'keepdiag':     keep the diags or not
% - 'flow':         Normalization: 1 for column-, or 2 for row-wise
% -------------------------------------------------------------------------
% OUTPUT
% - PDC:            PDC values: channel x  channel x freq x time x subj
% - Direction_Stats: a structure array with the following fields, calculated based on bootstrapping of directionalities: 
%                    LName: Layer name
%                    Median: Median of the bootstrap of directionalites
%                    CI:     confidence interval of the bootstraps
%                    PostStim:  stats calculated by comparing pre and post bootstrap distributions
%--------------------------------------------------------------------------
% Author: Elham Barzegaran, 09/2019
%
%% Default values

opt = ParseArgs(varargin, ...
        'nboot'         ,100, ...       % number of bootstraps
        'bootsize'      ,[],...         % size of the bootstrap draw: maximum = number of trials
        'ModOrds'       ,15, ...        % model order
        'ff'            ,.99, ...       % percent of variance too keep in STOK algorithm
        'measure'       ,'PDCnn', ...   % options: 'PDCnn' , 'PDC', 'sPDC', for more info see function: dynet_ar2pdc 
        'keepdiag'      ,1, ...         % 
        'flow'          ,2, ...          % 1 col, 2 row-wise normalization 
        'StatSide'      ,'right'...
        );

[~,nch,nt]       =    size(epochs{1}); % get the size of data
nsubj            =    numel(epochs);  % number of animals

if nch == 12% If it is recorded from two rois
    rnum    =   2;
else
    rnum    =   1;
end
lnum                =   round(nch/rnum); % number of layers/channels

%% (1) Estimate PDCs for each animal
disp('Calculating PDCs')

for S = 1:nsubj %     
    this_epochs = epochs {S}; % select the animal data

    %  estimate MVAR coeffs using stok
    M1 = max(max(abs(mean(this_epochs(:,:,tsec>0)))));
    M2 = max(max(abs(mean(this_epochs(:,:,tsec<0)))));
    M1/M2
        
     if M1/M2>10 % if the SNR is too high STOK cannot tune the C parameter properly
         fac = (0:3:30)*M2;
     else
        fac = 0;
     end
    
    parfor f = 1:numel(fac)
        Data              =    this_epochs(:,:,:)+rand(size(this_epochs))*fac(f);
        KF(f)             =    dynet_SSM_STOK(Data,opt.ModOrds,opt.ff);  
        Err(f)            =    sqrt(sum((KF(f).c(100:end)-KF(f).cu(100:end)).^2));
        
        M1 = max(max(abs(mean(Data(:,:,tsec>0)))));
        M2 = max(max(abs(mean(Data(:,:,tsec<0)))));
        SNR(f) = M1/M2;
    end
    
    % what SNR level to choose?
    if numel(fac)>1
        snr = min(find(diff(Err)==0));
        SNR(snr)
    else
        snr = 1;
    end
    
    %  estimate PDC based on MVAR coeffs
    PDC(:,:,:,:,S) =    dynet_ar2pdc(KF(snr),srate,fvec,opt.measure,opt.keepdiag,opt.flow); % We cannot keep the PDC bootstraps because of memory issues :(, check what are the other possibilities

end

%% (2) permute the PDC values/ for permutation test

parfor S = 1:nsubj
    % (3) calculate directionalities (Upwards) for each animal and permutation using PDC values

    for roi =1:rnum 
         lInd        =      (1:lnum)+(roi-1)*lnum;
         Data        =      PDC(lInd,lInd,:,:,S); 

         % Calculate directionality of the layers: 2 to 5
         Dir_layer   =    arrayfun(@(x) squeeze(mean(Data(1:x-1,x,:,:),1) - mean(Data(x+1:end,x,:,:),1)),2:lnum-1,'uni',false);
         Direct_layer(:,:,:,S,roi)  =    permute(cat(3,Dir_layer{:}),[3 1:2]);

         % Directionality calculated using full connectivity matrix/ Not layer-specific
         SD         =       size(Data);
         DataR      =       reshape(Data,SD(1),SD(2),prod(SD(3:4)));
         Dir_all    =       arrayfun(@(k) (sum(sum(triu(DataR(:,:,k)),1),2)-sum(sum(tril(DataR(:,:,k)),1),2))./(lnum*(lnum-1)/2),1:prod(SD(3:4)));
         Direct_all(:,:,S,roi)    =   reshape(Dir_all,SD(3:4));
    end
end

DL = permute(Direct_layer,[1 5 2:4]); DL = reshape(DL,size(DL,1)*size(DL,2),size(DL,3),size(DL,4),size(DL,5));
DA = permute(Direct_all,[4 1:3]);
 
clear DataR Data Direct_layer Dir_layer Dir_all Direct_all epochs;% clean up the memory
 
%% (4) draw the bootstraps and : bootstrap over animals and trials
disp('Bootstraping...')
if isempty(opt.bootsize) % Size of bootstrap draw
        opt.bootsize    =   nsubj;
end

if opt.nboot ==  1 % indicate the bootstrap draws
    Boots          =    1:nsubj;
else
    Boots          =    randi(nsubj,opt.nboot,opt.bootsize);
end

PDCs        =   size(PDC);
PDC_Mboot   =   zeros(PDCs(1:4));

for bs = 1:opt.nboot % bootstrap over subjects/animals
    if mod(bs,1)==100,     disp(['Bootsrap# ' num2str(bs)]);    end
    
    DL_boot(:,:,:,bs)   =   mean(DL(:,:,:,Boots(bs,:)),4);
    DA_boot(:,:,:,bs)   =   mean(DA(:,:,:,Boots(bs,:)),4);
    PDC_Mboot(:,:,:,:)  =   PDC_Mboot(:,:,:,:)+mean(PDC(:,:,:,:,Boots(bs,:)),5);
end

PDC_Mboot               =   PDC_Mboot/opt.nboot;
PDC(:,:,:,:,end+1)      =   PDC_Mboot;
PDC_Results.PDC = PDC;
PDC_Results.C = CB;
%% (5) Estimate frequency-wise pre-stimulus histogram and calculate stats based on that


% Direction_layers ((12+2) x fvec x tsec>0):
D_boot              =   cat(1,DL_boot,DA_boot);
Layers_Names = [];
for roi = 1:rnum
    Layers_Names    =   [Layers_Names cellfun(@(x) [x '_' ROIs{roi}],labels(2:5),'uni',false)];
end
for roi = 1:rnum
    Layers_Names    =   [Layers_Names ['All_' ROIs{roi}]];
end


parfor x = 1:size(D_boot,1)
    Direction_Stats(x)          =   prestim_bootstats(squeeze(D_boot(x,:,:,:)),tsec,.05/2,lower(opt.StatSide));% make the suprathreshold an input param
end

for x = 1:size(D_boot,1)
    Direction_Stats(x).LName    =   Layers_Names{x};
end


%% (6) Calculate cluster-level statistics

% will be the sum of sstats of the connected components 

%% (7) indicate significant clusters using cluster-stat of the random distribution

% make a histogram of Cluster stats of random perumted data and compare it
% to unpermuted data

end

