function [I,ind,indi,indj] = pairwiseMutualInformation(X)
X = X(:,:);
n = size(X,2);
I = zeros(n);
ind = triu(ones(n),1);
[indi,indj] = ind2sub([n n],find(ind(:)));
n = length(indi);
for k=1:n
    I(indi(k), indj(k)) = mutualInformation(X(:,indi(k)),X(:,indj(k)));
    I(indj(k), indi(k)) = I(indi(k), indj(k));
end
