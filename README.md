# pepkalc

## Protein electrostatics
Protonation is a ubiquitous and important process in biology. Protein folding, ligand recognition, enzyme catalysis, membrane potentials, and the energetics of cells depend on ionization and proton transfer. Charge-charge interactions are of special importance for Intrinsically Disordered Proteins, which are known to contain abnormally high numbers of consecutive charged amino acids. Consequently, a great theoretical effort has been devoted to elucidation of protein electrostatic interactions in unfolded state.

## The problem
Polypeptide sequence length is the single dominant factor hampering the effectiveness of currently available software tools for _de novo_ calculation of amino acid-specific protonation constants in disordered polypeptides.

## Our solution
**pepKalc** is robust simulation software for the comprehensive evaluation of protein electrostatics in unfolded state. Our software completely removes the limitations of the previously described Monte-Carlo approaches in the computation of protein electrostatics, by using a hybrid approach that effectively combines exact and mean-field calculations to rapidly obtain accurate results. Paired with a modern architecture CPUs, **pepKalc** is capable of evaluating protonation behavior for an arbitrary-size polypeptide in a sub-second time regime.

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
The Gaussian-chain parameter `A`. Default value `5.0`.

#### --gcb
The Gaussian-chain parameter `B`. Default value `7.5`.

#### --cutoff
The cutoff size for explicit interaction energy calculations. Default value `2`.  

#### --ncycles
The number of calculation super-cycles. Default value `3`.

#### --nooutput
Disable titration curve output.
