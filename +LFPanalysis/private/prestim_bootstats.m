function Stats = prestim_bootstats(Data,tsec,suprath,varargin)

% INPUTs:
    % Data: Is a freq x time x bootstraps matrix
    % tsec: Indicates the time (zero for stim onset) in ms
    % suprath: The supra-threshold for prestim significance indication

%%
opt = ParseArgs(varargin, ...
            'side',     'right' ... % 'right' or 'left' or 'both'
            );

if ~exist('suprath','var') || isempty(suprath)
    suprath = .01;
end


%% Estimate 95% confidence intervals and robust estimator of median

alpha   =   .05;

for F = 1:size(Data,1)
    for T = 1:size(Data,2)
        CI1(F,T)        =   quantile(Data(F,T,:),alpha/2);
        CI2(F,T)        =   quantile(Data(F,T,:),1-alpha/2);
        MedianEst(F,T)  =   median(Data(F,T,(Data(F,T,:)>=CI1(F,T))&(Data(F,T,:)<=CI2(F,T))));% robust estimator based on median of 95% CI
    end
end

Stats.CI        =   cat(3,CI1,CI2);
Stats.Median    =   MedianEst;


%%
presT           =   tsec<0;         % pre-stimulus time points for noise distribution estimation
presT(1:150)    =   false;          % allow for AR adjustment
posT            =   find(tsec>=0);  % post-stimulus time points

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
        M2      =       MedianEst(F,posT(pt));
        switch opt.side
            case 'right'
                h(F,pt)     =       double(M2>Th); % hyp test
                tstat(F,pt) =       sum(H>=Th & H<=M2)./numel(H); % significance level
                p(F,pt)     =       sum(H>M2)./numel(H); % p-value
                
            case 'left'           
                h(F,pt)     =       double(M2<Th);
                tstat(F,pt) =       sum(H<=Th & H>=M2)./numel(H);
                p(F,pt)     =       sum(H<M2)./numel(H);
                
            case 'both'
                h(F,pt)     =       double(M2>Th(1) || M2<Th(2));
                tstat(F,pt) =       (sum(H>=Th(1) & H<=M2)+sum(H<=Th(2) & H>=M2))./numel(H);
                p(F,pt)     =       (sum(H>=M2 & H>=TH(1))+sum(H<=M2 & H<=TH(2)))./numel(H); %%%?? HOW SHUOLD THIS BE?
                
        end    
        
    end
end

% post stim stats
Stats.PostStim.p        =   p;
Stats.PostStim.sstat    =   tstat;
Stats.PostStim.h        =   h;


end