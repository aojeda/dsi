# Distributed EEG source imaging toolbox for MATLAB/EEGLAB

![splash](https://github.com/aojeda/dsi/blob/master/doc/splash.png)

The DSI toolbox for is a metapackage that combines the [headModel](https://github.com/aojeda/headModel#headmodel-toolbox-for-matlabeeglab) and the [PEB+](https://github.com/aojeda/PEB) toolboxes into one coherent plugin for forward and inverse distributed source imaging within the [EEGLAB](https://sccn.ucsd.edu/eeglab/) environment.

## Install
1. [Download](https://github.com/aojeda/dsi/archive/master.zip)
2. Uncompress
3. Rename the folder as "EEGBrowser", and place it in your eeglab/plugins folder.
4. Run `eeglab` or `eeglab redraw` (if you don't want to clean your current EEG structure)

## Three functions to rule them all
* Forward problem solver: [pop_forwardProblem](https://github.com/aojeda/dsi/wiki/Forward-problem-solver)
* Inverse problem solver: [pop_inverseProblem](https://github.com/aojeda/dsi/wiki/Inverse-problem-solver)
* Visualization: [pop_eegbrowserx](https://github.com/aojeda/dsi/wiki/EEGBrowserX)
