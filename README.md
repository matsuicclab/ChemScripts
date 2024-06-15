# ChemScripts

scripts for input and output in computational chemical calculation

# setup
```
source <DIR>/bashrc # replace the text <DIR> to this directory according to your environment
```

# to use
## bash
```
$ atomnum2symb 'Ca'
20
$ smiles2xyz 'CCO'
$
```

## python
```
from pyg16.fchk import FCHK
fchk = FCHK('test.fchk')
fchk.calcElectronDensity([0, 0, 0])
```

# TODO
- g16log2value
    - eliminate title section
    - geometry
    - theory
- csvsmiles2xyz
- name2smiles
- bohr2ang
- ang2bohr




