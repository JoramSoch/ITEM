# ITEM

<h3>Inverse Transformed Encoding Models</h3>

This repository contains SPM-compatible MATLAB code for estimating inverse transformed encoding models (ITEM) based on first-level general linear models (GLMs) for functional magnetic resonance imaging (fMRI) data. The ITEM approach allows trial-wise linear decoding of discrete experimental conditions (classification) or continuous modulator variables (reconstruction) from multivariate fMRI signals.

An ITEM analysis would usually proceed in two steps:
1. Construct trial-wise design matrix <b>AND</b> estimate trial-wise response amplitudes using `ITEM_est_1st_lvl`. <br>
2. a) Classify discrete experimental conditions from trial-wise parameter estimates using `ITEM_dec_class` <b>OR</b> <br>
   b) Reconstruct continuous modulator variable from trial-wise parameter estimates using `ITEM_dec_recon`.

Type `help ITEM_fct_name` for information on input parameters of these functions. You may also use the review function by typing `ITEM_review` in order to check intermediate results at any time. The code in this repository references some functionality from the MACS toolbox [4]. If this toolbox is on the path, functions not starting by `ITEM_` are not required.

This software is in beta testing. Future improvements will include searchlight analysis, improved visualization, extended comments and a software manual. Preliminary documentation can be found in an abstract submitted to OHBM 2019 [1] as well as a preprint uploaded to <i>bioRxiv</i> [2] and currently under review at <i>NeuroImage</i>. In case of questions on the methodology or issues with the toolbox, please contact the corresponding author [3].

[1] https://ww5.aievolution.com/hbm1901/index.cfm?do=abs.viewAbs&subView=1&abs=3967 <br>
[2] https://www.biorxiv.org/content/10.1101/610626v1 <br>
[3] mailto:joram.soch@bccn-berlin.de?subject=ITEM%20toolbox%20inquiry <br>
[4] https://github.com/JoramSoch/MACS <br>
