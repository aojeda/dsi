function [I,ind,indi,indj] = pairwiseMutualInformation(X,normalize, kneig)
if nargin < 3, kneig=4;end
X = X(:,:);
n = size(X,2);
I = zeros(n);
ind = triu(ones(n),1);
[indi,indj] = ind2sub([n n],find(ind(:)));
nv = length(indi);
Iv = zeros(nv,1);
parfor k=1:nv
    Iv(k) = mutualInformation(X(:,indi(k)),X(:,indj(k)),normalize, kneig);
    %Iv(k) = mutualInformation_deprecated(X(:,indi(k)),X(:,indj(k)));
end
for k=1:nv
    I(indi(k), indj(k)) = Iv(k);
    I(indj(k), indi(k)) = Iv(k);
end