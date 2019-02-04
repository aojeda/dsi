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
addpath(fullfile(p,'hadModel'));
addpath(fullfile(p,'hadModel','resources'));
addpath(fullfile(p,'hadModel','+vis'));
addpath(fullfile(p,'hadModel','+geometricTools'));
addpath(genpath(fullfile(p,'hadModel','dependency')));
addpath(fullfile(p,'PEB'));
addpath(fullfile(p,'PEB','invChol'));
addpath(fullfile(p,'PEB','resources'));

h = uimenu( findobj(gcf, 'tag', 'tools'), 'label', 'Distributed spurce imaging');
hFp = uimenu( h, 'label','Forward problem');
hIp = uimenu( h, 'label','Inverse problem');

uimenu( hFp, 'label', 'Compute BEM forward model','callback','EEG = pop_forwardModel(EEG);');
uimenu( hFp, 'label', 'View head model','callback','hm=headModel.loadFromFile(EEG.etc.src.hmfile);hm.plot;');
uimenu( hFp, 'label', 'Documentation','callback','web(''https://github.com/aojeda/headModel'')');

uimenu( hIp, 'label', 'PEB+ source estimation','callback','EEG = pop_inverseSolution(EEG);');
uimenu( hIp, 'label', 'Source browser','callback','pop_eegbrowserx(EEG);');
uimenu( hIp, 'label', 'Move ROI source estimates to EEG.data','callback','try,EEG = moveSource2DataField(EEG);[ALLEEG EEG CURRENTSET]=eeg_store(ALLEEG, EEG);eeglab redraw;catch e, errordlg(e.message);end');
