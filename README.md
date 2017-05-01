# METSig
**(missing transerverse energy significance)**

Development of the MET significance variable for the 8 TeV dataset at CMS (CERN).  The variable is documented in Chapter 8 of [this paper](http://iopscience.iop.org/article/10.1088/1748-0221/10/02/P02006/meta;jsessionid=CBA6B45E51EF392FFC6D0A0F89B0A1D3.c4.iopscience.cld.iop.org).

## Features:
- Analyzes event information to estimate the MET covariance matrix.
- Algorithm is designed using Monte Carlo simulation and optimized on real data.  Optimization is conducted using a maximum likelihood fit.
- Outputs performance plots and optimal parameters for MET significance algorithm.

## Usage instructions:
- Ntuples from `METSigTuning/MakeNtuple` are indexed in `file_catalog.txt`.  Create using the script `./generate_catalog.sh`.
- To run likelihood fit for JER tuning and generate plots:
```
./DoFit <options>
```
- To generate plots and ROC curves without tuning jet energy resolutions:
```
./EvalSig <options>
```

### 76X Branches:
- 76XMINIAOD: For tuning on the Run2 dataset.  Currently set up for DY->dimuon events.

### 53X Branches:
- master: vanilla significance.
- fft_significance: contains FFT code as well as job parallelization scripts.  A merge of yimin_comment and nathan_fullshape.
- yimin_comment: Yimin Wang's final FFT code.
- systematics: for systematics evaluations.  Needs to be cleaned up, be cautious when using.
- nathan_pileup: generates ROC plots in bins of N pileup vertices
