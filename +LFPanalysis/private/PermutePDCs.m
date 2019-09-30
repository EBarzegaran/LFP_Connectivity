function PDC_perm = PermutePDCs(PDC)
% THIS function permutes PDC values in a PDC matrix, by permuting all the
% adjacency values, but keeping the diagonals constant
%%
    LNum = size(PDC,1);
    PDCs = size(PDC);
    RI = randperm(LNum*(LNum-1),LNum*(LNum-1));
    OI = 1:30; %This is in case of 6 Layers, should be updated later 
    OIorig = reshape(1:LNum*LNum,LNum,LNum); 
    RIorig = OIorig;
    ind = 1;
    for i = 1:LNum
        for j = 1:LNum
            if i~=j
                IndR(j,i) = RI(ind);
                IndO(j,i) = OI(ind);
                ind = ind+1;       
            end
        end
    end
    for i = 1:numel(IndR)
        RIorig(find(IndO==i)) = OIorig(find(IndR==i));
    end

    Datat = reshape(PDC,[LNum*LNum PDCs(3:end)]);
    PDC_perm = reshape(Datat(reshape(RIorig,LNum*LNum,1),:,:),PDCs); 
end