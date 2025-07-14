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
# generate smiles and structure
$ name='ethanol'
$ name2smiles "$name" | smiles2xyz -O test.xyz

# calc energy with Gaussian
$ xyz2gjf -O test_0.gjf --link0 '%chk=test_0.chk' --route '# opt b3lyp/6-31G(d)' --charge 0 test.xyz
$ g16 test_0.gjf
$ g16log2value --SCF-energy test_0.log
-155.033799257

# calc energy of cation with Gaussian
$ xyz2gjf -O test_1.gjf --link0 '%chk=test_1.chk' --route '# opt b3lyp/6-31G(d)' --charge 1 test.xyz
$ g16 test_1.gjf
$ g16log2value --SCF-energy test_1.log
-154.664421473

# calc DeltaE in kcal/mol
$ (g16log2value --SCF-energy test_1.log ; g16log2value --SCF-energy test_0.log) | energy2energy --from a.u. --to kcal/mol --format '%.10f' | tr '\n' ' ' | awk '{print $1-$2}'
231.788
```

## python
```
from pyg16.fchk import FCHK
fchk = FCHK('test.fchk')
fchk.calcElectronDensity([0, 0, 0])

from pyg16.cube import visualizeCubeSlice
visualizeCubeSlice(cubeFile='test.cub', outFile='test_slice.pdf')
```

# Are you in trouble?
## Error1
```
ChemScripts/bash/abstract/data2data: line 15: exec: {ChemScript_data2data_stdout_fd}: not found
```
This is probably because the version of bash you are using is old.
Check your bash version with /bin/bash --version.
bash4.1 or higher is recommended.


# TODO
- g16log2value
    - eliminate title section
    - geometry
    - theory
- g16log2xyz (wrapper of g16log2value --opt-geometry)
- g16log2csv
- csvsmiles2xyz
- name2xyz (For complexes that are difficult to express in SMILES notation)
- glog2prep (extend of mkresp)
- xyz2prep
- mksolv
- bkill-all -> bkill 0
- jdel-all
- sdf2mol
- pdb2pdbs
- xyz2gms
- mkmdin
- mdcrd2crd (ambpdb)
- xyz2xyz (consideration on isotope)
- file2seqfiles
- files2file
- void2file

- error handling for pdbid2pdb, name2sdf
- Cube class (generate cube with node data of original cube)
- convert CRLF to LF before processing


