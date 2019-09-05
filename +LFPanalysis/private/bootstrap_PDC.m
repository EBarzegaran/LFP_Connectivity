function  [PDC, C] = bootstrap_PDC(epochs,nboot,ModOrds,ff,srate,fvec,measure,keepdiag,flow)

% INPUTS:
    %   epochs: trial x channel x time
    %   nboot: number of bootstraps
    %   ModOrds: model order of STOK
%%

[ntrl,nch,nt] = size(epochs);

%% draw the bootstraps
nboot = 1;
Boots = 1:ntrl;%randi(ntrl,nboot,round(ntrl/2));

% (1) Build the bootstrap distribution
for b = 1:nboot
    if mod(b,10)==1, disp(['Bootsrap# ' num2str(b)]); end
    %  estimate MVAR coeffs using stok
    KF = dynet_SSM_STOK(epochs(Boots,:,:),ModOrds,ff);  
    %  estimate PDC based on MVAR coeffs
    CB(:,b) = KF.c;
    PDC(:,:,:,:,b) = dynet_ar2pdc(KF,srate,fvec,measure,keepdiag,flow);
end


C = mean(CB,2);
%check the C parameter to see how it fits



% (2) estimate frequency-wise prestim histogram
% arrayfun
% histcount

end