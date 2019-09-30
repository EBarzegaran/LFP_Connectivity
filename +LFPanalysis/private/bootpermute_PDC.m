function  [PDC_Results, Direction_Stats]  =   bootpermute_PDC  (epochs,  srate,  fvec,   tsec, labels, ROIs,  varargin)

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
% - labels:     a cell array of labels of the layers
% - ROIs:       if it is a multisite recording, the name of ROIs
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
        'StatSide'      ,'right',...
        'DoPerm'        ,false, ...
        'NPerm'         ,[], ...
        'animID'        ,[] ...
        );

[~,nch,nt]       =    size(epochs{1}); % get the size of data
nsubj            =    numel(epochs);  % number of animals

if nch == 12% If it is recorded from two rois
    rnum    =   2;
else
    rnum    =   1;
end
lnum                =   round(nch/rnum); % number of layers/channels

if opt.DoPerm % does either bootstrapping or permutation
    if isempty(opt.NPerm)
        opt.NPerm = 500;
    end
    opt.nboot = 1;
else
    opt.NPerm = 1;
end

% Layer Names
Layers_Names = [];
for roi = 1:rnum
    Layers_Names    =   [Layers_Names cellfun(@(x) [x '_' ROIs{roi}],labels(2:5),'uni',false)];
end
for roi = 1:rnum
    Layers_Names    =   [Layers_Names ['All_' ROIs{roi}]];
end
    
%% (1) Estimate PDCs for each animal
disp('Calculating PDCs')

for S = 1:nsubj %     
    this_epochs = epochs {S}; % select the animal data
    disp([num2str(round(S/nsubj*100)) '%']);
    %  estimate MVAR coeffs using stok
    M1 = max(max(abs(mean(this_epochs(:,:,tsec>0)))));
    M2 = max(max(abs(mean(this_epochs(:,:,tsec<0)))));
    SNR_orig(S) = M1/M2;
        
     if M1/M2>10 % if the SNR is too high STOK cannot tune the C parameter properly
         fac = (0:3:30)*M2;
     else
        fac = 0;
     end
    
    parfor f = 1:numel(fac)
        Data              =    this_epochs(:,:,:)+rand(size(this_epochs))*fac(f);
        KF(f)             =    dynet_SSM_STOK(Data,opt.ModOrds,opt.ff);  
        Err(S,f)          =    sqrt(sum((KF(f).c(100:end)-KF(f).cu(100:end)).^2));

        
        M1 = max(max(abs(mean(Data(:,:,tsec>0)))));
        M2 = max(max(abs(mean(Data(:,:,tsec<0)))));
        SNR(S,f) = M1/M2;
    end
    
    % what SNR level to choose?
    if numel(fac)>1
        snr = min(find(diff(squeeze(Err(S,:)))==0));
        snr_idx(S) = snr;
    else
        snr = 1;
    end
    
    CB(S,:)           =    KF(snr).c;% C parameter for later visualization
    
    %  estimate PDC based on MVAR coeffs
    PDC(:,:,:,:,S) =    dynet_ar2pdc(KF(snr),srate,fvec,opt.measure,opt.keepdiag,opt.flow); % We cannot keep the PDC bootstraps because of memory issues :(, check what are the other possibilities

end
%% SNR data visualization
if false
    FIG = figure;
    hold on;
    ylim([-50 800]);
    for S = 1:nsubj
        SP(S) = plot(10*log10(SNR(S,:)),Err(S,:),'linewidth',1.5,'marker','+');
        set(gca, 'XDir','reverse');
        line([SNR(S,snr_idx(S)) SNR(S,snr_idx(S))],[-50 800],'color',get(SP(S),'color'),'linestyle','--','linewidth',1.5);
    end
    
    xlabel('10log10(SNR)');
    ylabel('(C-C_{adjusted})^2');
    if ~isempty(opt.animID)
        legend(SP,opt.animID);
    end
    set(gca,'fontsize',12);
    set(FIG,'unit','inch','position',[5 5 10 5],'color','w');
    export_fig(FIG,'SNR_adjustment','-pdf'); close ;
end

%% (2) permute the PDC values/ for permutation test
disp('Permutation test...')
NPerm       =   opt.NPerm;
for perm    =   1:NPerm
    % generate labels for permutation
    if mod(perm,10)==0
        disp([num2str(round(perm/NPerm*100)) '%']);
    end
    %tic
    for S = 1:nsubj
        % (3) calculate directionalities (Upwards) for each animal and permutation using PDC values
        
        % permute the same for each subject according to the generated permutations
        for roi =1:rnum 
            lInd        =      (1:lnum)+(roi-1)*lnum; % which layers to select, according to ROI
            
            % Do not permute for the first sample, otherwise permute
            if perm==1
                PDCperm  =      PDC(lInd,lInd,:,:,S);
            else
                PDCperm  =      PermutePDCs(PDC(lInd,lInd,:,:,S));
            end
             
             Data        =      PDCperm; 
             %Data        =      PDC(lInd,lInd,:,:,S); 

             % Calculate directionality of the layers: 2 to 5
             Dir_layer   =    arrayfun(@(x) squeeze(mean(Data(1:x-1,x,:,:),1) - mean(Data(x+1:end,x,:,:),1)),2:lnum-1,'uni',false);
             Direct_layer(:,:,:,S,roi,perm)  =    permute(cat(3,Dir_layer{:}),[3 1:2]);

             % Directionality calculated using full connectivity matrix/ Not layer-specific
             SD         =       size(Data);
             DataR      =       reshape(Data,SD(1),SD(2),prod(SD(3:4)));
             Dir_all    =       arrayfun(@(k) (sum(sum(triu(DataR(:,:,k)),1),2)-sum(sum(tril(DataR(:,:,k)),1),2))./(lnum*(lnum-1)/2),1:prod(SD(3:4)));
             Direct_all(:,:,S,roi,perm)    =   reshape(Dir_all,SD(3:4));
        end
    end
    %toc
end
DL = permute(Direct_layer,[1 5 2:4 6]); DL = reshape(DL,size(DL,1)*size(DL,2),size(DL,3),size(DL,4),size(DL,5),size(DL,6));
DA = permute(Direct_all,[4 1:3 5]);
 
clear DataR Data Direct_layer Dir_layer Dir_all Direct_all epochs;% clean up the memory
 
%% (4) draw the bootstraps and : bootstrap over animals and trials
if ~opt.DoPerm
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

    % (5) Estimate frequency-wise pre-stimulus histogram and calculate stats based on that

    % Direction_layers ((8+2) x fvec x tsec>0):
    D_boot                  =   cat(1,DL_boot,DA_boot);
    
    parfor x = 1:size(D_boot,1)
        Direction_Stats(x)          =   prestim_bootstats(squeeze(D_boot(x,:,:,:)),tsec,.05/2,lower(opt.StatSide));% make the suprathreshold an input param
    end

    for x = 1:size(D_boot,1)
        Direction_Stats(x).LName    =   Layers_Names{x};
    end
else
    
    %% Permutation test
    % prepare output
    PDC(:,:,:,:,end+1)      =   mean(PDC,5);
    PDC_Results.PDC         =   PDC;
    PDC_Results.C           =   CB;
    
    % (4) Permutation test with cluster-level stats
    D              =   cat(1,DL,DA);
    parfor x = 1:size(D,1)
        Direction_Stats(x)          =   permute_stats(D(x,:,:,:),tsec,.05/2,lower(opt.StatSide));
    end
    
    
end


end

