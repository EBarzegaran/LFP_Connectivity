function  [PDC, C]  =   bootstrap_PDC  (epochs,  srate,  fvec,   tsec,  varargin)

% This function calculates PDC based on STOK algorithm and using
% bootstrapping it estimnate the significant results.

% INPUTS:
    %   epochs: trial x channel x time
    %   nboot: number of bootstraps
    %   ModOrds: model order of STOK

% OUTPUT:
    %   PDC:
    %   C:
    %   Direction:
    
%     
%% default values

opt = ParseArgs(varargin, ...
        'nboot'         ,100, ...       % number of bootstraps
        'bootsize'      ,[],...         % size of the bootstrap draw: maximum = number of trials
        'ModOrds'       ,15, ...        % model order
        'ff'            ,.99, ...       % percent of variance too keep in STOK algorithm
        'measure'       ,'PDCnn', ...   % options: 'PDCnn' , 'PDC', 'sPDC', for more info see function: dynet_ar2pdc 
        'keepdiag'      ,1, ...         % 
        'flow'          ,2 ...          % 1 col, 2 row-wise normalization 
        );

[ntrl,nch,nt]       =    size(epochs); % get the size of data

if isempty(opt.bootsize) % size of bootstrap draw
    opt.bootsize    =   round(ntrl/3);
end

%% (1) draw the bootstraps and calculate PDCs

if opt.nboot ==  1
    Boots          =    1:ntrl;
else
    Boots          =    randi(ntrl,opt.nboot,opt.bootsize);
end

% Build the bootstrap distribution of PDC estimates
for b = 1:opt.nboot
    if mod(b,1)==0,     disp(['Bootsrap# ' num2str(b)]);    end
    
    %  estimate MVAR coeffs using stok
    KF             =    dynet_SSM_STOK(epochs(Boots(b,:),:,:),opt.ModOrds,opt.ff);  
    
    %  estimate PDC based on MVAR coeffs
    CB(:,b)        =    KF.c;
    PDC(:,:,:,:,b) =    dynet_ar2pdc(KF,srate,fvec,opt.measure,opt.keepdiag,opt.flow);
    
end
C   =    mean(CB,2); % C parameter indicates the dynamics of PDC over time

%% (2) calculate directionality (Upwards) distributions using PDC dists

if nch == 12% If it is recorded from two rois
    rnum = 2;
else
    rnum = 1;
end

lnum = round(nch/rnum); % number of layers/channels

 for roi =1:rnum 
     lInd        =      (1:lnum)+(roi-1)*lnum;
     Data        =      PDC(lInd,lInd,:,:,:); 
     
     % Calculate directionality of the layers
     Dir_layer          =    arrayfun(@(x) squeeze(mean(Data(1:x-1,x,:,:,:),1) - mean(Data(x+1:end,x,:,:,:),1)),1:lnum,'uni',false);
     Direct_layer{roi}  =    permute(cat(4,Dir_layer{:}),[4 1:3]);
     
     % Directionality calculated using full connectivity matrix/ Not layer-specific
     SD         =       size(Data);
     DataR      =       reshape(Data,SD(1),SD(2),prod(SD(3:5)));
     Dir_all    =       arrayfun(@(k) (sum(sum(triu(DataR(:,:,k)),1),2)-sum(sum(tril(DataR(:,:,k)),1),2))./(lnum*(lnum-1)/2),1:prod(SD(3:5)));
     Direct_all{roi}    =   reshape(Dir_all,SD(3:5));

 end

 DL = cat(1,Direct_layer{:});
 DA = permute(cat(4,Direct_all{:}),[4 1:3]);
 
 clear DataR Data Direct_layer Dir_all Direct_all;% clean up the memory
%% (3) estimate frequency-wise pre-stimulus histogram

% PDC for each channel pair (12 x 12 x fvec x tsec>0)

% Direct_layers (12 x fvec x tsec>0): Refered as DL

% Direct_all (2 x fvec x tsec>0): Refered as DA

prestim_bootstats(squeeze(DA(1,:,:,:)),tsec)

end

