# ChemScripts

scripts for input and output in computational chemical calculation

# setup
1. download by `git clone` or `wget`
```
git clone https://github.com/matsuicclab/ChemScripts
or
wget https://github.com/matsuicclab/ChemScripts/archive/master.tar.gz
tar xpvf master.tar.gz
```

2. set environment variables (append to ~/.bashrc)
```
source <DIR>/bashrc # replace the text <DIR> to this directory according to your environment
```

# to use
## bash
```
$ atomnum2symb 'Ca'
20
$ smiles2xyz 'CCO'
$ g16log2value --SCF-energy test.log
-76.4341246819
```

## python
```
from pyg16.fchk import FCHK
fchk = FCHK('test.fchk')
fchk.calcElectronDensity([0, 0, 0])

from pyg16.cube import visualizeCubeSlice
visualizeCubeSlice(cubeFile='test.cub', outFile='test_slice.pdf')
```

# TODO
- g16log2value
    - eliminate title section
    - geometry
    - theory
- g16log2xyz (wrapper of g16log2value --opt-geometry)
- g16log2csv
- xyz2gjf
- csvsmiles2xyz
- name2smiles
- name2xyz (For complexes that are difficult to express in SMILES notation)
- pdbid2pdb
- glog2prep (extend of mkresp)
- smiles2prep
- bkill-all -> bkill 0
- jdel-all




