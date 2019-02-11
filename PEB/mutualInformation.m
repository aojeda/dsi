function I = mutualInformation(x,y,normalize,k)
% Use the JIDT library to compute Kraskov's estimate of mutual information 
% between continuous random variables x and y.
%
% For this function to work it should be located next to infodynamics.jar, 
% otherwise download the jar file from: http://jlizier.github.io/jidt/
%
%% Example:
%%  Two variables that share very little information:
%   x=randn(100000,1);
%   y=randn(size(x));
%   I = mutualInformation(x,y);
%   r = corr(x,y);
%   disp([I r])
%
%%  Two variables that share information in a complicated manner:
%   x=randn(100000,1);
%   y=x.^2;
%   I = mutualInformation(x,y);
%   r = corr(x,y);
%   disp([I r])
%
%%  Two variables that are linearly related:
%   x=randn(100000,1);
%   y=x+0.1*randn(100000,1);
%   I = mutualInformation(x,y);
%   r = corr(x,y);
%   disp([I r])
%
% If you use this function for a publication please cite as:
% Joseph T. Lizier, "JIDT: An information-theoretic toolkit for studying 
% the dynamics of complex systems", Frontiers in Robotics and AI 1:11, 2014;
% doi:10.3389/frobt.2014.00011 (pre-print: arXiv:1408.3270)
%
% Author of the wrapper: Alejandro Ojeda, Neural Engineering and Translation
%                        Labs, University of California San Diego, 2019

persistent miCalc
if isempty(miCalc)
    p = fileparts(which('mutualInformation.m'));
    javaaddpath(fullfile(p, 'infodynamics.jar'));
    miCalc=javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');
    miCalc.initialise(1,1); % univariate calculation
end
if nargin < 3
    normalize = false;    % Set to true to zscore input signals
end
if nargin < 4
    k=4;                    % Nearest neighbours for KSG estimator
end
if normalize
    normalize = 'true';
else
    normalize = 'false';
end
miCalc.setProperty('k', num2str(k)); 
miCalc.setProperty('NORMALISE', normalize);
miCalc.setObservations(x(:), y(:));
I = miCalc.computeAverageLocalOfObservations();
end
