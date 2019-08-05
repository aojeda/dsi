function [H, Delta, blocks, indG, indV] = buildAugmentedLeadField(hm)
% Augments the lead field matrix with artifact scalp maps projections.
% 
% Inputs:
%   hm: is a head model object
%   M: is a dictionary of artifact scalp map projections (each column is a scalp map)
%   chanlocs: is the channel locations structure that corresponds to each column of M 
%
% Outputs:
%   H: [L A] is the lead field matrix augmented with a dictionary of artifact scalp maps
%   indG are the indices of the brain sources
%   indV are the indices of the artifact sources

Nroi = length(hm.atlas.label);
[Ny,Ng] = size(hm.K);
c = load('Artifact_dictionary.mat');
templateFile = fullfile(fileparts(which('headModel.m')),'resources',[c.templateName '.mat']);
template = headModel.loadFromFile(templateFile);
%--
Nic = size(c.A,2);
A = zeros(Ny,Nic);
for ic=1:Nic
    F = scatteredInterpolant(template.channelSpace,c.A(:,ic),'linear','linear');
    A(:,ic) = F(hm.channelSpace);
end
A = [A eye(Ny)];
% A = eye(Ny);
%--
Nv = size(A,2);
norm_K = norm(hm.K);
L = hm.K/norm_K;
Delta = hm.L/norm_K;
H = [L A];
H = bsxfun(@rdivide,H,sqrt(sum(H.^2)));
blocks = hm.indices4Structure(hm.atlas.label);
if size(hm.K,2) == 3*size(hm.cortex.vertices,1)
    H(:,(1:Ng)) = -H(:,(1:Ng));
    Delta = kron(eye(3),Delta);
    blocks = logical(kron(eye(3),blocks));
    blocks = [[blocks;false(Nv,Nroi*3)] [false(Ng,Nv);diag(true(Nv,1))]];
else
    blocks = [[blocks;false(Nv,Nroi)] [false(Ng,Nv);diag(true(Nv,1))]];
end
Delta = blkdiag(Delta,eye(Nv));
indG = (1:Ng)';
indV = Ng+(1:Nv)';
end