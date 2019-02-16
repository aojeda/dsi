function EEG = pop_srcConnectivity(EEG, time, connType, normalize, neig)
if nargin < 2
    time = EEG.times([1 end]);
end
if nargin < 3
    connType = 'te';
end
if nargin < 4
    normalize = false;
end
if nargin < 5
    neig = 4;
end
ind = EEG.times > time(1) & EEG.times < time(2);
[Nroi,~, Ntrials] = size(EEG.etc.src.act(:,ind,:));
C = zeros([Nroi Nroi, Ntrials]);

X = permute(EEG.etc.src.act(:,ind,:),[2 1 3]);
if strcmp(connType,'te')
    theta = zeros([Nroi Nroi, Ntrials]);
    for k=1:Ntrials
        disp(['Trial ' num2str(k) '/' num2str(Ntrials)]);
        [C(:,:,k), theta(:,:,k)] = conditionalTransferEntropy(X(:,:,k), -1, normalize, neig);
    end
    EEG.etc.conn = struct('measure','conditional transfer entropy','C',C,'theta',theta);
elseif strcmp(connType,'mi')
    [C(:,:,1),ind,indI,indJ] = pairwiseMutualInformation(X(:,:,1));
    ind = find(ind(:)==1);
    Cv = zeros(length(ind),Ntrials);
    tmp = C(:,:,1);
    Cv(:,1) = tmp(ind);
    for k=2:Ntrials
        disp(['Trial ' num2str(k) '/' num2str(Ntrials)]);
        C(:,:,k) = pairwiseMutualInformation(X(:,:,k));
        tmp = C(:,:,k);
        Cv(:,k) = tmp(ind);
    end
    EEG.etc.conn = struct('measure','mutual information','C',C,'Cv',Cv,'indI',indI,'indJ',indJ);
end
disp('done!')
    