# ERAASR

ERAASR (Estimation and Removal of Array Artifacts via Sequential principal components Regression) is an algorithm for removing electrical stimulation artifacts from multielectrode array recordings. 

Note: if you are looking for the version that works with Neuropixels recordings, please contact me (Github handle at Stanford). This version should be published shortly and differs meaningfully from the version in this repo, although the concept is similar.

Journal of Neural Engineering (accepted manuscript): http://iopscience.iop.org/article/10.1088/1741-2552/aaa365

Preprint on bioRxiv: https://www.biorxiv.org/content/early/2017/09/07/185850. 

Source code (Matlab): https://github.com/djoshea/eraasr/

## Installation:
- If you want to use the example dataset, first install [Git LFS (large file storage)](https://github.com/git-lfs/git-lfs/wiki/Installation) which will pull down the `exampleDataTensor.mat` when you clone the repo. Alternatively you can download the [`exampleDataTensor.mat`](https://github.com/djoshea/eraasr/raw/master/exampleDataTensor.mat) file directly.
- Clone the repo: `git clone https://github.com/djoshea/eraasr.git`
- Add `eraasr` directory to your Matlab path. Only the root path is necessary, e.g. `addpath '/path/to/eraasr'`, since all of the functions are inside the `ERAASR` package. [Matlab packages](https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) operate like namespaces.

## Instructions:

Walk through the code in [`example.m`](https://github.com/djoshea/eraasr/blob/master/example.m) to get started with the example dataset. You can copy this script and change parameters as needed for your own datasets. If you leave the `showFigures` parameter to true, ERAASR will generate a variety of diagnostic figures that can show you how each stage of the algorithm is working on your data. When you are satisfied, setting `showFigures` to false will speed up the code considerably. This example script also demonstrates how to reinsert the cleaned data back into each trial, high-pass filter the cleaned data, and threshold and extract spiking waveforms using tools provided with ERAASR.
