function val = resampleROI(value,vertices, faces, roi, npoints)
f = faces(all(ismember(faces,roi),2),:);
convexHull = [[...
    vertices(f(:,1),:)+vertices(f(:,2),:);...
    vertices(f(:,1),:)+vertices(f(:,3),:);...
    vertices(f(:,2),:)+vertices(f(:,3),:)]/2;...
    (vertices(f(:,1),:)+vertices(f(:,2),:)+vertices(f(:,3),:))/3];
points = [vertices(roi,:);convexHull];
if npoints > size(points,1)
    sample = randi(size(points,1),npoints,1);
else
    sample = randperm(size(points,1),npoints);
end
X1 = points(sample,:);
X2 = vertices(roi,:);
n1 = size(X1,1);
n2 = size(X2,1);
X1 = reshape(repmat(X1,n2,1),[n1 n2 3]);
X2 = permute(reshape(repmat(X2,n1,1),[n2 n1 3]),[2 1 3]);
W = 1./(eps+sqrt(sum((X1-X2).^2,3)));
W = bsxfun(@rdivide, W, eps+sum(W,2));
val = W*value;
end
