% eegplugin_dsi() - Distributed EEG source imaging toolbox for EEGLAB.
% Usage:
%   >> eegplugin_dsi();
%
% Author: Alejandro Ojeda, NEATLabs, University of California San Diego, 2019
%
% See also: eeglab()

function vers = eegplugin_dsi(fig,try_strings, catch_strings)
vers = 'DIS1.0.0';
p = fileparts(which('eegplugin_dsi'));
addpath(fullfile(p,'headModel'));
addpath(fullfile(p,'headModel','resources'));
addpath(genpath(fullfile(p,'headModel','dependency')));
addpath(fullfile(p,'RSBL'));
addpath(fullfile(p,'RSBL','invChol'));
addpath(fullfile(p,'RSBL','resources'));
addpath(fullfile(p,'Connectivity'));

h = uimenu( findobj(gcf, 'tag', 'tools'), 'label', 'Distributed source imaging');
hFp = uimenu( h, 'label','Forward problem');
hIp = uimenu( h, 'label','Inverse problem');
uimenu( h, 'label','Create my artifact dictionary','callback','CreateArtifactDictionary');

uimenu( hFp, 'label', 'Compute BEM forward model','callback','EEG = pop_forwardModel(EEG);');
uimenu( hFp, 'label', 'View head model','callback','hm=headModel.loadFromFile(EEG.etc.src.hmfile);hm.plot;');
uimenu( hFp, 'label', 'Documentation','callback','web(''https://github.com/aojeda/headModel'')');

uimenu( hIp, 'label', 'RSBL source estimation','callback','EEG = pop_rsbl(EEG);');
uimenu( hIp, 'label', 'Source browser','callback','pop_eegbrowserx(EEG);');
uimenu( hIp, 'label', 'Move ROI source estimates to EEG.data','callback','try,EEG = moveSource2DataField(EEG);[ALLEEG EEG CURRENTSET]=eeg_store(ALLEEG, EEG);eeglab redraw;catch e, errordlg(e.message);end');
