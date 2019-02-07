function I = mi(x,y)
[px, xi]   = ksdensity(x);
[py, yi]   = ksdensity(y);
[pxy, xyi] = ksdensity([x(:),y(:)]);

% Evaluate pdfs in the common support
n = 100;
cs = unique([xi yi]);
cs = linspace(min(cs),max(cs),n);
fx = griddedInterpolant(xi,px);
fy = griddedInterpolant(yi,py);
px = fx(cs);px(px<eps)=0;
py = fy(cs);py(py<eps)=0;

X = reshape(xyi(:,1),[30 30]);
Y = reshape(xyi(:,2),[30 30]);
pxy = reshape(pxy,[30 30]);
[Sx,Sy] = meshgrid(cs);
pxy = interp2(X,Y,pxy,Sx,Sy);
pxy(isnan(pxy)) = 0;

% Compute MI as: sum_x sum_y Pxy*(log(Pxy)-log(Px)-log(Py))
logPxy = log(pxy);
logPxy(isinf(logPxy)) = prctile(logPxy(~isinf(logPxy)),1);
logPxPy = bsxfun(@plus,log(px'),log(py));
logPxPy(isinf(logPxPy)) = prctile(logPxPy(~isinf(logPxPy)),1);

I = pxy.*(logPxy - logPxPy);
I(isnan(I)|isinf(I)) = 0;
dx = cs(2)-cs(1);
dy = cs(2)-cs(1);
I = sum(sum(I*dy)*dx);
