function [T, theta] = conditionalTransferEntropy(X, theta, normalize, k)
% Use the JIDT library to compute Kraskov's estimate of transfer entropy
% between continuous random variables x and y, conditioned on a set of
% variables.
%
% Here we approximate the delays at which the TE is calculated by finding
% the maximum of the autocorrelation function between pair of variables.
%
% For this function to work it should be located next to infodynamics.jar,
% otherwise download the jar file from: http://jlizier.github.io/jidt/
%
%% Examples:
%%  x(k-10) -> y (k)
%   n = 1000; 
%   lag = 10;
%   x=randn(n,1);
%   y=circshift(x,lag)+0.1*randn(n,1);
%   z = 0.5*x+0.5*randn(n,1);
%   T = conditionalTransferEntropy([x,y,z]);
%   r = corr([x,y,z]);
%   subplot(121);imagesc(T);colorbar;
%   subplot(122);imagesc(r);colorbar;
%
%%  z(k-10) -> y(k), z(k-10) -> x(k)
%   n = 1000;
%   lag = 10;
%   z=randn(n,1);
%   y=0.6*circshift(z,lag)+0.1*randn(n,1);
%   x=0.4*circshift(z,lag)+0.1*randn(n,1);
%   T = conditionalTransferEntropy([x,y,z]);
%   r = corr([x,y,z]);
%   subplot(121);imagesc(T);colorbar;
%   subplot(122);imagesc(r);colorbar;
%

%%  z(k-10) -> y(k)^2 <-x(k-5)^3
%   n = 1000;
%   lag = 10;
%   z=randn(n,1);
%   x=randn(n,1);
%   y=0.6*circshift(z,lag).^2+0.6*circshift(x,5).^2+0.1*randn(n,1);
%   [T, theta] = conditionalTransferEntropy([x,y,z]);
%   r = corr([x,y,z]);
%   subplot(121);imagesc(T);colorbar;
%   subplot(122);imagesc(r);colorbar;
%%
%%
% If you use this function for a publication please cite as:
% Joseph T. Lizier, "JIDT: An information-theoretic toolkit for studying
% the dynamics of complex systems", Frontiers in Robotics and AI 1:11, 2014;
% doi:10.3389/frobt.2014.00011 (pre-print: arXiv:1408.3270)
%
% Author of the wrapper: Alejandro Ojeda, Neural Engineering and Translation
%                        Labs, University of California San Diego, 2019
%%
persistent isInJavaPath
if isempty(isInJavaPath)
    clear -java;
    p = fileparts(which('transferEntropy.m'));
    javaaddpath(fullfile(p, 'infodynamics.jar'));
    isInJavaPath = true;
end
if nargin < 2
    theta = -1;             % Number of samples the source is ahead of the target.
                            % If theta=-1 we use the mutual information function
                            % to find it.
end
if nargin < 3
    normalize = false;      % If true we zscore all input signals
end
if nargin < 4
    k=4;                    % Nearest neighbours for KSG estimator
end
if normalize
    normalize = 'true';
else
    normalize = 'false';
end
[n,m] = size(X);
teCalc=javaObject('infodynamics.measures.continuous.kraskov.ConditionalTransferEntropyCalculatorKraskov');
teCalc.setProperty('k', num2str(k));
teCalc.setProperty('NORMALISE', normalize);
T = zeros(m);
if theta == -1
    theta = zeros(m);
    findLag = true;
else
    findLag = false;
end
I = zeros(n,1);
ind = triu(true(m),1);
[indI, indJ] = ind2sub([m,m],find(ind(:)));

Nc = length(indI);
prc = round(Nc*(0.1:0.1:1));
for c=1:Nc
    if findLag
        for k=1:n
            I(k) = mutualInformation(X(:,indJ(c)),circshift(X(:,indI(c)),k-1));
        end
        
        % Find a significant delay
        z = zscore(I);
        d = find(z > norminv(1-(0.005/n)));
        if ~isempty(d)
            [~,loc] = max(z(d));
            d = d(loc)-1;
            theta(indJ(c),indI(c)) = d;
            if theta(indJ(c),indI(c)) > n/2
                d = (n-theta(indJ(c),indI(c)));
                tmp = indI(c);
                indI(c) = indJ(c);
                indJ(c) = tmp;
                theta(indJ(c),indI(c)) = d;
                theta(indI(c),indJ(c)) = 0;
            end
        else
            continue
        end
    end
    teCalc.initialise(...
        1,1, ...                                        % Destination embedding length (Schreiber k=1) and delays
        1,1, ...                                        % Source embedding length (Schreiber l=1) and delays
        abs(theta(indJ(c),indI(c))), ...                % Source-destination delay of 1 (default)
        ones(1,m-2), ...                                % Embedding lengths for each conditional variable
        ones(1,m-2)*abs(theta(indJ(c),indI(c))), ...    % Embedding delays for each conditional variable
        ones(1,m-2)*abs(theta(indJ(c),indI(c))) ...     % Conditional-destination delays for each conditional variable
        );
    teCalc.setObservations(X(:,indI(c)), X(:,indJ(c)), X(:,setdiff(1:m,[indI(c) indJ(c)])));
    T(indJ(c),indI(c)) = teCalc.computeAverageLocalOfObservations();
    
    % Progress indicatior
    prc_c = find(prc==c);
    if ~mod(c,10)
        fprintf('.');
    elseif ~isempty(prc_c)
        fprintf('%i%%',(prc_c)*10)
    end
end
clear -java
end
