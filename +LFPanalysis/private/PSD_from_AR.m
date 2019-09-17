function PSD = PSD_from_AR(AR,p,srate,freqs)

% AR is channel x channel x p x time
% p is model order
%%
nodes = size(AR,1);
freqs = freqs(:);
%--------------------------------------------------------------------------
    Z        = repmat(exp(-2*pi*1i*freqs/srate),[1 p]).^repmat((1:p),[numel(freqs) 1]);
    A        = repmat(eye(nodes), [1 1 length(freqs) size(AR,4)]);

    for k = 1:p
        tmp  = repmat(-AR(:, :, k,:), [1 1 length(freqs) 1]);
        A    = A + bsxfun(@times,tmp,reshape(Z(:,k),1,1,[]));
    end

    NCov = 1;

    PSD = (1/srate*NCov)./abs(repmat(eye(size(A,1)),1,1,size(A,3),size(A,4))+A).^2;
    
        
end