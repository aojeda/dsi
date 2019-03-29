function EEG = pop_srcConnectivity(EEG, time, connType, normalize, neig, useFull)
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
if nargin < 6
    useFull = false;
end
ind = EEG.times > time(1) & EEG.times < time(2);
[Nroi,~, Ntrials] = size(EEG.etc.src.act(:,ind,:));
C = zeros([Nroi Nroi, Ntrials]);
hm = headModel.loadFromFile(EEG.etc.src.hmfile);

if strcmp(connType,'te')
    X = permute(EEG.etc.src.act(:,ind,:),[2 1 3]);
    theta = zeros([Nroi Nroi, Ntrials]);
    for k=1:Ntrials
        disp(['Trial ' num2str(k) '/' num2str(Ntrials)]);
        [C(:,:,k), theta(:,:,k)] = conditionalTransferEntropy(X(:,:,k), -1, normalize, neig);
    end
    EEG.etc.conn = struct('measure','conditional transfer entropy','C',C,'theta',theta);
elseif strcmp(connType,'mi')
    if useFull
        X = permute(EEG.etc.src.actFull(EEG.etc.src.indG,ind,:),[2 1 3]);
        if normalize
            X = bsxfun(@rdivide, X, sqrt(sum(sum(X.^2))));
        end
        [C(:,:,1),ind,indI,indJ] = pairwiseMutualInformationROI(X(:,:,1), hm);
        ind = find(ind(:)==1);
        Cv = zeros(length(ind),Ntrials);
        tmp = C(:,:,1);
        Cv(:,1) = tmp(ind);
        for k=2:Ntrials
            disp(['Trial ' num2str(k) '/' num2str(Ntrials)]);
            C(:,:,k) = pairwiseMutualInformationROI(X(:,:,k), hm);
            tmp = C(:,:,k);
            Cv(:,k) = tmp(ind);
        end
    else
        X = EEG.etc.src.act(:,ind,:);
        if normalize
            X = bsxfun(@rdivide, X, sqrt(sum(sum(X.^2))));
        end
        C = pairwiseMutualInformation(X(:,:)', normalize);
        X_surrogate = X(:,:)';
        ntnt = size(X_surrogate,1);
        for k=1:size(X_surrogate,2)
            X_surrogate(:,k) = X_surrogate(randperm(ntnt,ntnt),k);
        end
        C_surrogate = pairwiseMutualInformation(X_surrogate, normalize);
        X = permute(X,[2 1 3]);
        disp(['Trial 1/' num2str(Ntrials)]);
        [Ck(:,:,1),ind,indI,indJ] = pairwiseMutualInformation(X(:,:,1),normalize);
        for k=2:Ntrials
            disp(['Trial ' num2str(k) '/' num2str(Ntrials)]);
            Ck(:,:,k) = pairwiseMutualInformation(X(:,:,k),normalize);
        end
    end
    EEG.etc.conn = struct('measure','mutual information','C',C,'C_surrogate',C_surrogate,'Ck',Ck,'ind',ind,'indI',indI,'indJ',indJ);
end
disp('done!')
    