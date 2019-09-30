function Stats = permute_stats(Data,tsec,suprath,side)
%
% This function calculasted the PDC/directionality statistics based on permutation test:
% random distribution

% Syntax: Stats = permute_stats(Data,tsec,suprath,side)
%-------------------------------------------------------------
%  INPUTs:
%  - Data: Is a freq x time x numberofperms matrix: first permutation
%          element should be the unpermuted data
%  - tsec: Indicates the time (zero for stim onset) in ms
%  - suprath: The supra-threshold for prestim significance indication
%  - side: direction of hypothesis testing: ['right'] or 'left' or 'both'  
%--------------------------------------------------------------------------
% Author: Elham Barzegaran, 10/2019
%
%% Default values
if ~exist('side','var') || isempty(suprath)% 'right' or 'left' or 'both'
    side = 'right';
end

if ~exist('suprath','var') || isempty(suprath)
    suprath = .01;
end
Data = squeeze(Data);
Nperm = size(Data,3);

%% For each time point and frequency, indicate if PDC/directionality is significant based on random distribution
for np = 1:Nperm
    for F = 1:size(Data,1)
        for T  = 1:size(Data,2)
            if strcmpi(side,'right')
                p(F,T) = sum(Data(F,T,:)>Data(F,T,np))./Nperm;
            elseif strcmpi(side,'left')
                p(F,T) = sum(Data(F,T,:)<Data(F,T,np))./Nperm;
            end
        end
    end
    % cluster level stats
    [CS,Clusters] = findclust(p,Data(:,:,1),suprath);
    CS_dist(np) = CS(1);
    
    if np==1
        Stats.p = p;
        Stats.Clusters = Clusters;
        Stats.ClustersSS = CS;
    end

    
end

% final cluster level p-value
Stats.ClustersP = arrayfun(@(x) sum(CS_dist>Stats.ClustersSS(x))/Nperm,1:numel(Stats.ClustersSS),'uni',false);

end
