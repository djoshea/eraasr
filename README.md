# ERAASR

ERAASR (Estimation and Removal of Array Artifacts via Sequential principal components Regression) is an algorithm for removing electrical stimulation artifacts from multielectrode array recordings. Read more in the preprint up on bioRxiv at https://www.biorxiv.org/content/early/2017/09/07/185850. 

Source code (Matlab): https://github.com/djoshea/eraasr/

Installation:
- If you want to use the example dataset, first install [Git LFS (large file storage)](https://github.com/git-lfs/git-lfs/wiki/Installation) which will pull down the `exampleDataTensor.mat` when you clone the repo. Alternatively you can download the [`exampleDataTensor.mat`](https://github.com/djoshea/eraasr/raw/master/exampleDataTensor.mat) file directly.
- Clone the repo: `git clone https://github.com/djoshea/eraasr.git`
- Add `eraasr` directory to your Matlab path. Only the root path is necessary, e.g. `addpath '/path/to/eraasr'`, since all of the functions are inside the `ERAASR` package. [Matlab packages](https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) operate like namespaces.
- Walk through the code in [`example.m`](https://github.com/djoshea/eraasr/blob/master/example.m) to get started
