function I = mutualInformation(x,y, pdfType)
% Compute the mutual information (in nats) between continuous random variables x and
% y.
%
% Example:
%  Two variables that share very little information:
%   x=randn(100000,1);
%   y=randn(size(x));
%   I = mutualInformation(x,y);
%   r = corr(x,y);
%   disp([I r])
%
%  Two variables that share information in a complicated manner:
%   x=randn(100000,1);
%   y=x.^2;
%   I = mutualInformation(x,y);
%   r = corr(x,y);
%   disp([I r])
%
%  Two variables that are linearly related:
%   x=randn(100000,1);
%   y=x+0.1*randn(100000,1);
%   I = mutualInformation(x,y)
%   r = corr(x,y);
%   disp([I r])
%
% Author: Alejandro Ojeda, Neural Engineering and Translation Labs, University of California San Diego, 2019

if nargin < 3, pdfType = 'ksdensity';end
if ~any(ismember({'ksdensity','hist'},lower(pdfType)))
    pdfType = 'ksdensity';
end

x = x/(eps+std(x(:)));
y = y/(eps+std(y(:)));
if strcmpi(pdfType,'ksdensity')
    [px, xi]   = ksdensity(x);
    [py, yi]   = ksdensity(y);
    [pxy, xyi] = ksdensity([x(:),y(:)]);
else
    [px, xi]   = hist(x,100);
    [py, yi]   = hist(y,100);
    [pxy, xyi] = hist3([x(:),y(:)],[30 30]);
    pxy = pxy(:);
    [Sx,Sy] = meshgrid(xyi{1},xyi{2});
    xyi = [Sx(:),Sy(:)];
end

% Evaluate pdfs in the common support
n = 100;
cs = unique([xi yi]);
cs = linspace(min(cs),max(cs),n);
fx = griddedInterpolant(xi,px);
fy = griddedInterpolant(yi,py);
px = fx(cs);px(px<eps)=0;
py = fy(cs);py(py<eps)=0;

nxy = sqrt(length(pxy));
pxy = reshape(pxy,[nxy,nxy]);
X = reshape(xyi(:,1),[nxy nxy]);
Y = reshape(xyi(:,2),[nxy nxy]);
[Sx,Sy] = meshgrid(cs);
pxy = interp2(X,Y,pxy,Sx,Sy,'linear',0);
pxy(pxy<eps) = 0;

% Normalize the pdfs
px = px/nansum(px);
py = py/nansum(py);
pxy = pxy/nansum(pxy);

% Compute MI as: sum_x sum_y Pxy*(log(Pxy)-log(Px)-log(Py))
logPxy = log(pxy);
logPxPy = bsxfun(@plus,log(px'),log(py));
logPxy(isinf(logPxy)) = min(logPxy(~isinf(logPxy)));
logPxPy(isinf(logPxPy)) = min(logPxPy(~isinf(logPxPy)));
I = pxy.*(logPxy - logPxPy);
I(isinf(I)) = nan;
dx = cs(2)-cs(1);
dy = cs(2)-cs(1);
I = nansum(I(:)*dy*dx);
if I<0 && strcmpi(pdfType,'hist')
    I = mutualInformation(x,y, 'ksdensity');
end