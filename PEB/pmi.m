function I = pmi(X)
X = X(:,:);
n = size(X,1);
I = zeros(n);
ind = triu(ones(n));
[indi,indj] = ind2sub([n n],find(ind(:)));
n = length(indi);
for k=1:n
    I(indi(k), indj(k)) = mi(X(indi(k),:),X(indj(k),:));
    I(indj(k), indi(k)) = I(indi(k), indj(k));
end
