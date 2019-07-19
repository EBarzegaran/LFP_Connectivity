function REV = LSK_REV2(data, M_ord, l)
%--------------------------------------------------------------------------
% [Mattia - 07.05.2016]
% Computation of Relative Error Covariance (REV):
% REV = MSE/MSY
% MSE: mean square of prediction error (residuals)
% MSY: variance of the signal (mean squared signal Y)
% 1-REV can be used as a measure of goodness-of-fit
%--------------------------------------------------------------------------
% [1]	Alois Schlögl (2000). The Electroencephalogram and the Adaptive
%       Autoregressive Model: Theory and Applications
%==========================================================================

% gp adapted for LSK, 05/2019

%KF = SSM_KF_STTK(data,M_ord,l);
KF = dynet_SSM_STOK(data,M_ord,l);

% % Option 1: Y observed and Y estimated with PAR coefficients
% [nTr, nChan, nSamp, ] = size(data);
% Yest = zeros(nTr, nChan, nSamp); 
% for p = 1:M_ord
%     for t = M_ord+1:nSamp
%         Yest(:,:,t) = Yest(:,:,t) + (KF.AR(:,:,p,t)*data(:,:,t-p)')';        %
%     end
% end
% Yest(:,:,1:M_ord) = [];
% data(:,:,1:M_ord) = [];
% KF.PY(:,:,1:M_ord) = [];
% 
% MSE = mean((data(:)-Yest(:)).^2);
% MSY = mean(data(:).^2);
% 
% MSE = var(data(:)-Yest(:));
% MSY = var(data(:));
% 
% REV = MSE / MSY;
% end option 1

% Option 2: Y observed and one step prediction error
data(:,:,1:M_ord) = [];
KF.PYe(:,:,1:M_ord) = [];
MSE2 = mean((KF.PY(:)).^2);
MSY2 = mean(data(:).^2);
REV = MSE2 / MSY2;

end    









