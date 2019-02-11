function TE = transferEntropy(x,y,lag,normalize,k)
% Use the JIDT library to compute Kraskov's estimate of transfer entropy 
% between continuous random variables x -> y.
%
% For this function to work it should be located next to infodynamics.jar, 
% otherwise download the jar file from: http://jlizier.github.io/jidt/
%
%% Examples:
%%  x(k-10) -> y (k)
%   lag = 10;
%   x=randn(1000,1);
%   y=circshift(x,lag)+0.1*randn(size(x));
%   TE = transferEntropy(x,y,lag);
%   r = corr(x,y);
%   disp([TE r])
%
%%  x(k-10) -> y(k)^2
%   lag = 10;
%   x=randn(1000,1);
%   y=circshift(x,lag).^2+0.1*randn(size(x));
%   TE = transferEntropy(x,y,lag);
%   r = corr(x,y);
%   disp([TE r])
%
%%  x is unrelated to y
%   lag = 10;
%   x=randn(1000,1);
%   y=randn(size(x));
%   TE = transferEntropy(x,y,lag);
%   r = corr(x,y);
%   disp([TE r])
%
%%
% If you use this function for a publication please cite as:
% Joseph T. Lizier, "JIDT: An information-theoretic toolkit for studying 
% the dynamics of complex systems", Frontiers in Robotics and AI 1:11, 2014;
% doi:10.3389/frobt.2014.00011 (pre-print: arXiv:1408.3270)
%
% Author of the wrapper: Alejandro Ojeda, Neural Engineering and Translation
%                        Labs, University of California San Diego, 2019
%%
persistent teCalc
if isempty(teCalc)
    p = fileparts(which('transferEntropy.m'));
    javaaddpath(fullfile(p, 'infodynamics.jar'));
    teCalc=javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');
    teCalc.initialise(1);
end
if nargin < 3
    lag = 1;                % Number of samples x is ahead of y
end
if nargin < 4
    normalize = false;      % Set to true to zscore input signals
end
if nargin < 5
    k=4;                    % Nearest neighbours for KSG estimator
end
if normalize
    normalize = 'true';
else
    normalize = 'false';
end
teCalc.setProperty('DELAY',num2str(lag));
teCalc.setProperty('k', num2str(k)); 
teCalc.setProperty('NORMALISE', normalize);
teCalc.setObservations(x(:), y(:));
TE =teCalc.computeAverageLocalOfObservations();
end
