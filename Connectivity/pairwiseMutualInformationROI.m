function [I,ind,indi,indj] = pairwiseMutualInformationROI(X, hm)
roi = hm.indices4Structure(hm.atlas.label);
X = X(:,:);
n = size(roi,2);
Xc = zeros(64,size(X,1),n);
for k=1:n
    Xc(:,:,k) = geometricTools.resampleROI(X(:,roi(:,k))',hm.cortex.vertices, hm.cortex.faces, find(roi(:,k)), 64);
end

I = zeros(n);
ind = triu(ones(n),1);
[indi,indj] = ind2sub([n n],find(ind(:)));
nv = length(indi);
Iv = zeros(n);
tic
parfor k=1:nv
%     ni = sum(roi(:,indi(k)));
%     nj = sum(roi(:,indj(k)));
%     if ni<nj
%         xi = geometricTools.resampleROI(X(:,roi(:,indi(k)))',hm.cortex.vertices, hm.cortex.faces, find(roi(:,indi(k))), nj);
%         xj = X(:,roi(:,indj(k)))';
%     else
%         xj = geometricTools.resampleROI(X(:,roi(:,indj(k)))',hm.cortex.vertices, hm.cortex.faces, find(roi(:,indj(k))), ni);
%         xi = X(:,roi(:,indi(k)))';
%     end
    xi = Xc(:,:,indi(k));
    xj = Xc(:,:,indj(k));
    Iv(k) = mutualInformation(xi(:),xj(:),true);
end
for k=1:nv
    I(indi(k), indj(k)) = Iv(k);
    I(indj(k), indi(k)) = Iv(k);
end
toc