function [I,ind,indi,indj] = pairwiseMutualInformation(X,normalize)
X = X(:,:);
n = size(X,2);
I = zeros(n);
ind = triu(ones(n),1);
[indi,indj] = ind2sub([n n],find(ind(:)));
nv = length(indi);
Iv = zeros(nv);
parfor k=1:nv
    Iv(k) = mutualInformation(X(:,indi(k)),X(:,indj(k)),normalize);
end
for k=1:nv
    I(indi(k), indj(k)) = Iv(k);
    I(indj(k), indi(k)) = Iv(k);
end