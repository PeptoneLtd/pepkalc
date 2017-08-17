import argparse
from scipy.optimize import curve_fit
from scipy.special import erfc
import numpy as np
global aadict

# Our amino acid conversion utilities
three = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
         "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
one = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
       "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
aadict = {}

for i in range(len(three)):
    aadict[three[i]] = one[i]
    aadict[one[i]] = three[i]

def one2three2one(aa):
    """
    Amino acid 1to3 conversion
    """
    return(aadict[aa.upper()])

def get_diagonal(matrix, n=0):
    z = np.zeros(matrix.shape, int)
    for i in range(len(z) - n):
        z[i, i + n] = True
    return matrix[z > 0]

def smallmatrixlimits(ires, cutoff, len):
    ileft = max(1, ires - cutoff)
    iright = min(ileft + 2 * cutoff, len)
    if iright == len:
        ileft = max(1, iright - 2 * cutoff)
    return (ileft, iright)

def smallmatrixpos(ires, cutoff, len):
    resi = cutoff + 1
    if ires < cutoff + 1:
        resi = ires
    if ires > len - cutoff:
        resi = min(len, 2 * cutoff + 1) - (len - ires)
    return resi

def hasselbalch(pH, pK, nH):
    return (10 ** (pK - pH)) ** nH / (1. + (10 ** (pK - pH)) ** nH)

def k(ion, T, eps):
    return np.sqrt((8.0 * np.pi * ion * 0.2003457539870666) / (eps * T * 0.0019872041))

def W1(r, ion, T, eps):
    kappa = k(ion, T, eps)
    x = kappa * r / np.sqrt(6.0)
    return 332.286 * np.sqrt(6.0 / np.pi) * (1 - np.sqrt(np.pi) * x * np.exp(x ** 2) * erfc(x)) / (eps * r)

def WP(ze, c, ion, T, eps):
    kappa = k(ion, T, eps)
    r = 11.8 / pow(c, 1.0 / 3.0)  # in Angstroms
    return ze * ze * np.exp(-kappa * r) / (eps * r)

def W(ds, ion, eps, weight=None):
    k = np.sqrt(ion) / 3.08
    ds[ds == 0] = 1
    tmp = 332.286 / (eps * ds) * np.exp(- k * ds)
    return np.average(tmp, axis=2, weights=weight)

def w2logp(x, R, T):
    return x * 4181.2 / (R * T * np.log(10))

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--sequence", default='nMDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEAc')
parser.add_argument("--temperature", default=283.15)
parser.add_argument("--ionicstrength", default=0.0)
parser.add_argument("--epsilon", default=83.83)
parser.add_argument("--gca", default=5)
parser.add_argument("--gcb", default=7.5)
parser.add_argument("--cutoff", default=2)
parser.add_argument("--ncycles", default=3)
args = parser.parse_args()

print 'pepKalc paramteres:'
print (args)

# Variables init
seq = args.sequence
temp = float(args.temperature)
ion = float(args.ionicstrength)
eps = float(args.epsilon)
gca = float(args.gca)
gcb = float(args.gcb)
cutoff = int(args.cutoff)
ncycles = int(args.ncycles)

pK0 = {"n": 7.5, "C": 8.6, "D": 4.0, "E": 4.4,
       "H": 6.6, "K": 10.4, "R": 12.0, "Y": 9.6, "c": 3.5}

pos = np.array([i for i in xrange(len(seq)) if seq[i] in pK0.keys()])
sites = ''.join([seq[i] for i in pos])

l = np.array([abs(pos - pos[i]) for i in range(len(pos))])
d = gca + np.sqrt(l) * gcb
tmp = W1(d, ion, temp, eps)

R = 8.314472
pHstep = 0.2
pHs = np.arange(1.0, 13.01, pHstep)

neg = np.array([i for i in range(len(pos)) if seq[pos[i]] in 'CDEYc'])

#print('distancematrix: ', ds.shape)
I = np.diag(np.ones(len(pos)))
tmp[I == 1] = 0
ww = w2logp(tmp, R, temp) / 2

chargesempty = np.zeros(pos.shape[0])
if len(neg):
    chargesempty[neg] = -1

pK0s = [pK0[c] for c in sites]
nH0s = [0.9 for c in sites]

titration = np.zeros((len(pos), len(pHs)))
smallN = min(2 * cutoff + 1, len(pos))
smallI = np.diag(np.ones(smallN))

alltuples = [[int(c) for c in np.binary_repr(i, smallN)]
             for i in xrange(2 ** (smallN))]
gmatrix = [np.zeros((smallN, smallN)) for p in range(len(pHs))]

if len(sites) <= 2 * cutoff + 1:
    ncycles = 1

titration = np.array([[hasselbalch(pHs[p], pK0s[i], nH0s[i])
                       for i in range(len(pos))] for p in range(len(pHs))]).transpose()

for icycle in range(ncycles):
    # print icycle + 1
    fractionhold = titration.transpose()

    for ires in range(1, len(pos) + 1):
        (ileft, iright) = smallmatrixlimits(ires, cutoff, len(pos))
        resi = smallmatrixpos(ires, cutoff, len(pos))

        for p in range(len(pHs)):
            fraction = fractionhold[p].copy()
            fraction[ileft - 1: iright] = 0
            charges = chargesempty + fraction
            ww0 = np.diag(np.dot(ww, charges) * 2)
            gmatrixfull = ww + ww0 + pHs[p] * I - np.diag(pK0s)
            gmatrix[p] = gmatrixfull[ileft - 1: iright, ileft - 1: iright]

        E_all = np.array([sum([10 ** -(gmatrix[p] * np.outer(c, c)).sum()
                               for c in alltuples]) for p in range(len(pHs))])
        E_sel = np.array([sum([10 ** -(gmatrix[p] * np.outer(c, c)).sum()
                               for c in alltuples if c[resi - 1] == 1]) for p in range(len(pHs))])
        titration[ires - 1] = E_sel / E_all

    sol = np.array([curve_fit(hasselbalch, pHs, titration[i], [
                   pK0s[i], nH0s[i]])[0] for i in range(len(pK0s))])
    (pKs, nHs) = sol.transpose()

first = 0
if not seq[0] == 'n':
    first = 1

# Print the header
print '%s, %s, %s' % ('RES','dpKa', 'n')
for i in range(len(sol)):
    print '%s, %f, %f' % (seq[pos[i]] + str(pos[i] + first), pKs[i] - pK0s[i], nHs[i]) #pKs[i]#, pKs[i] - pK0s[i], nHs[i]
    outfile = open(seq[pos[i]] + str(pos[i] + first) + '_titration.dat', 'w')
    for p in range(len(pHs)):
        outfile.write('%7.7f, %7.7f\n' % (pHs[p], titration[i][p]))
    outfile.close()
