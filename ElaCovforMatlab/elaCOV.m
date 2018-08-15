function [R,S] = elaCOV(X,Cr, Cv)

% INPUT
% X: nxp observation matrix 
% Cr: penlaty parameter for the estimation of correlation coefficient
%      matrix.
% Cv: penlaty parameter for the estimation of variances 
%
%     For large dimension Cr = 2 and Cv = 2 are suggested.
%
% OUTPUT
% R: estimated correlation coefficient matrix of pxp
% S: estimated covariance matrix of pxp
%
% Remarks: if you find the Matlab code is not fast enough, an C++ code is
%          available as explained below.


[n,p]=size(X);
x = X - repmat(mean(X),n,1);

nCr=length(Cr);

[S,S0]= elaCore(x,Cr+1.0e-10);  % if you find this is still slow please use the next sentence 

% please compile myMatrixCore.cpp first 
% [S,S0]= myMatrixCore(x,Cr+1.0e-10); 

R = reshape(S,p,p,nCr);

V = var(X);
Vk = zeros(p,1);
W = std(X.^2)/sqrt(n);
for i = 1:p
    Ii = (V<=V(i)+Cv*W(i))&(V>=V(i)-Cv*W(i));
    Vk(i) = mean(V(Ii));
end
V = sqrt(Vk);
S = V*V'.*R;

end

