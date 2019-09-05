function [M, Y, C, L] = buildArtifactDictionary(dataFolders, jProgressBar)
srate = 128;
template = headModel.loadDefault;
M = [];

%% Select files
files = [];
for k=1:length(dataFolders)
    tmp = pickfiles(dataFolders{k},'.set');
    files = char(files, tmp);
end
files(1,:) = [];
n = size(files,1);

for subject=1:n
    file = deblank(files(subject,:));
    [filePath,fileName,ex] = fileparts(file);
    if exist(fullfile(filePath,[fileName '_Colin27',ex]),'file')
        continue;
    end
    EEG = pop_loadset(file);
    
    %% Run ICA
    EEG.data = double(EEG.data);
    EEG = pop_rejcont(EEG, 'elecrange',1:EEG.nbchan , 'freqlimit',[20 40] ,...
        'threshold',10,'epochlength',0.5,'contiguous',4,...
        'addlength',0.25,'taper','hamming');
    EEG = pop_reref( EEG, []);
    EEG = pop_resample( EEG, srate);
    EEG = pop_eegfiltnew(EEG, [],0.5,846,1,[],0);  % 0.5 Hz highpass
    EEG = pop_eegfiltnew(EEG, 48,52,212,1,[],0);
    EEG = pop_eegfiltnew(EEG, 58,62,212,1,[],0);
        
    EEG = pop_runica(EEG, 'extended',1,'icatype','binica');
    f = pickfiles(pwd,'bias_after_adjust'); for d=1:size(f,1), delete(deblank(f(d,:)));end
    f = pickfiles(pwd,'.sc'); for d=1:size(f,1), delete(deblank(f(d,:)));end
    f = pickfiles(pwd,'.sph');for d=1:size(f,1), delete(deblank(f(d,:)));end
    f = pickfiles(pwd,'.wts');for d=1:size(f,1), delete(deblank(f(d,:)));end
    
    %% Co-register with template
    EEG = pop_forwardModel(EEG);
    hm = headModel.loadFromFile(EEG.etc.src.hmfile);
    delete(EEG.etc.src.hmfile);
    
    %% Interpolate ICs in the channel space of the template
    for ic = 1:size(EEG.icawinv,2)
        F = scatteredInterpolant(hm.channelSpace,EEG.icawinv(:,ic));
        M = [M F(template.channelSpace)];     %#ok
        % Plot IC on the skin of the template to make sure that the interpolation went well 
        % template.plotOnModel(randn(size(template.cortex.vertices),1),M(:,end)); 
    end
    jProgressBar.setValue(100*subject/n);
end

%% T-SNE
Y = tsne(M','Distance','correlation','NumDimensions',2,'Verbose',1,'Options',struct('MaxIter',500,'OutputFcn',[],'TolFun',1e-10),'Perplexity',50);

%% K_Means
nc = 8;
[L, C] = kmeans(M',nc,'distance','sqeuclidean', 'Replicates',9);
C = C';
end