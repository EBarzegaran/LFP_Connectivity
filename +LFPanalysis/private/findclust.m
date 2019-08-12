function [CS,Clusters] = findclust(PData,Data,Suprath)
    % find clusters and return the indices and their summary statistics
    
    BData = PData<Suprath;% binary: 1 for suprathreshold significant and 0 nonsignificant
    Data(isnan(Data))=0;
    % concat 0 to first and last elements of array
    BSize = size(BData);
    CC = bwconncomp(BData);
    Clusters = CC.PixelIdxList;
    CS = cellfun(@(x) sum(Data(x)),Clusters);
    [CS,Ind] = sort(CS,'descend');
    Clusters = Clusters(Ind);
    
%     %%    % Build the adjacency matrix for connected component extraction
%     Inds = reshape(1:numel(BData),size(BData));
% 
%     A = zeros(max(Inds(:))); 
%     for r = 1:BSize(1)
%         for c = 1:BSize(2)
%             if r>1, A(Inds(r,c),Inds(r-1,c))=1; end
%             if r<BSize(1), A(Inds(r,c),Inds(r+1,c))=1; end
%             
%             if c>1, A(Inds(r,c),Inds(r,c-1))=1; end
%             if c<BSize(2), A(Inds(r,c),Inds(r,c+1))=1; end
%             
%         end
%     end
%     A = max(A,A');
%     A(1:length(A)+1:end)=1;
%     A2 = A;
%     A2(BData==0,:)=0;
%     A2(:,BData==0)=0;
%     %% find cluster statistics
%     Clust = conncomp(graph(A2));
%     for C = 1:max(Clust)
%         CS(C)= sum(Data(Clust==C));
%     end
%     
%     [CS IndC]= sort(CS,'descend');
%     for C = 1:max(Clust)
%         Clust2(Clust==IndC(C))=C;
%     end
%     Clust2 = reshape(Clust2,BSize).*BData;
end