# pepkalc

## Installation
Clone this repository,
```
https://github.com/ktamiola/pepkalc.git
```
and install Python dependencies (you will need `pip` utility for that),
```
pip install scipy numpy
```

## Usage
Call **pepkalc** like any other Python script, parsing command-line paramters, e.g.
```
python pepkalc.py --sequence DDD
```
will perform pKa and Hill parameter estimations for `DDD` polypeptide.

## Parameters
**pepkalc** accepts the following input parameters:

#### --sequence
One-letter amino acid sequence following FASTA convention. Please use `n` and `c` to include N- and C-Terminus in your calculations. Default value `nMDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEAc`.

#### --temperature
The temperature in `K`. Default value `283.15`.

#### --ionicstrength
The ionic strength in `M`. Default value `0.0`.

#### --epsilon
The dielectric permeability of solvent. Default value `83.83` (assuming aqueous solution).

#### --gca
The Gaussian-chain parameter `A`. Default value `7.5`.

#### --gcb
The Gaussian-chain parameter `B`. Default value `5.0`.

#### --cutoff
The cutoff size for explicit interaction energy calculations. Default value `2`.  

#### --ncycles
The number of calculation super-cycles. Default value `3`.
