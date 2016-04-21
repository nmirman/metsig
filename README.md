# Analysis code for MET Significance.

# Usage
   * To compile, run make.
   * To read ntuples from METSigTuning/MakeNtuple, use ./generate_catalog.sh to create a text file with all ntuple file names.
   * To run likelihood fit for JER tuning and generate plots, run ./DoFit.
   * To generate plots and ROC curves without tuning JERs, run ./EvalSig.

# 76X Branches:
   * 76XMINIAOD: For tuning on the Run2 dataset.  Currently set up for DY->dimuon events.

# 53X Branches:
   * master: vanilla significance.
   * fft_significance: should have FFT code as well as job parallelization scripts.  A merge of yimin_comment and nathan_fullshape.
   * yimin_comment: Yimin's final FFT code.
   * systematics: has code for systematics evaluations.  Needs to be cleaned up, be cautious when using.
   * nathan_pileup: for making ROC plots in bins of N vertices
