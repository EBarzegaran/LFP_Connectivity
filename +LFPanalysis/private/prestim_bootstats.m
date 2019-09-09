function Stats = prestim_bootstats(Data,tsec,suprath,varargin)

% INPUTs:
    % Data: is a freq x time x bootstraps matrix
    % tsec: indicates the time (zero for stim onset) in ms

%%
opt = ParseArgs(varargin, ...
            'side',     'right' ... % 'right' or 'left' or 'both'
            );

if ~exist('suprath','var') || isempty(suprath)
    suprath = .01;
end
tic

%%
presT = tsec<0;         % pre-stimulus time points for noise distribution estimation
presT(1:150) = false;   % allow for AR adjustment

posT = find(tsec>=0);

for F = 1:(size(Data,1))
    H = Data(F,presT,:);
    H       =   sort(H(:));
    
    switch opt.side
        case 'right'            
            Th      =   H(round((1-suprath)*length(H)));
            
        case 'left'           
            Th      =   H(round(suprath*length(H)));
            
        case 'both'
            Th(1)   =   H(round((1-suprath/2)*length(H)));
            Th(2)   =   H(round((suprath/2)*length(H)));
        otherwise
            error ('uncorrect side option...')
    end
    
    for pt = 1:numel(posT)
        H2      =       sort(squeeze(Data(F,posT(pt),:)));
        switch opt.side
            case 'right'
                h(F,pt)     =       double(mean(H2)>Th);
                tstat(F,pt) =       sum(H>=Th & H<=mean(H2))./numel(H);
                p(F,pt)     =       sum(H>mean(H2))./numel(H);
                
            case 'left'           
                h(F,pt)     =       double(mean(H2)<Th);
                tstat(F,pt) =       sum(H<=Th & H>=mean(H2))./numel(H);
                p(F,pt)     =       sum(H<mean(H2))./numel(H);
                
            case 'both'
                h(F,pt)     =       double(mean(H2)>Th(1) || mean(H2)<Th(2));
                tstat(F,pt) =       (sum(H>=Th(1) & H<=mean(H2))+sum(H<=Th(2) & H>=mean(H2)))./numel(H);
                p(F,pt)     =       (sum(H>=mean(H2) & H>=TH(1))+sum(H<=mean(H2)& H<=TH(2)))./numel(H); %%%?? HOW SHUOLD THIS BE?
                
        end    
        
    end
end



%% Cluster-based correction?


end