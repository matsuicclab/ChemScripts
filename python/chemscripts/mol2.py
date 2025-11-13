import re

from chemscripts.molecule import Molecule

class Mol2:
    def __init__(self, filePath):
        mol2Data = open(filePath, mode='r').readlines()
        keyIdxs = [i for i,line in enumerate(mol2Data) if '@<TRIPOS>' in line]
        dataList = [mol2Data[i:f] for i,f in zip(keyIdxs, keyIdxs[1:]+[len(mol2Data)])]
        dataDict = {re.sub('@<TRIPOS>', '', data[0]).strip():data[1:] for data in dataList}
        self.__mol2Data = mol2Data
        self.__dataDict = dataDict
        self.__numAtom = int(dataDict['MOLECULE'][1].split()[0])
        self.__numBond = int(dataDict['MOLECULE'][1].split()[1])

        xyzBlock = [line.split()[1:5] for line in dataDict['ATOM'][:self.__numAtom]]
        xyzBlock = '\n'.join(['{} {} {} {}'.format(re.sub('[0-9]+','',s),x,y,z) for s,x,y,z in xyzBlock])
        xyzBlock = '{}\nloaded from {} by chemscripts\n{}'.format(self.__numAtom,filePath,xyzBlock)
        self.__molecule = Molecule(xyzBlock=xyzBlock)

    def giveXYZBlock(self, unit='Angstrom', elementSymbol=True, comment='', atomfilter=None, xyzformat=None):
        return self.__molecule.giveXYZBlock(unit=unit, elementSymbol=elementSymbol, comment=comment, atomfilter=atomfilter, xyzformat=xyzformat)

    def modifyBondTypeForAntechamber(self, octetcharge=0):
        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        from rdkit.Chem import ResonanceMolSupplier, ResonanceFlags

        def generateModifiedBondData(xyzBlock, octetcharge=0):
            mol = Chem.MolFromXYZBlock(xyzBlock)
            rdDetermineBonds.DetermineBonds(mol, charge=octetcharge)
            suppl = ResonanceMolSupplier(mol, ResonanceFlags.KEKULE_ALL)
            bondTypeList = [[b.GetBondType() for b in m.GetBonds()] for m in suppl]
            bondTypeList = np.array(bondTypeList)
            # bondTypeが各極限構造式で変わる結合(n%1!=0)については'ar'に置換
            bondType = [str(int(n)) if n%1==0 else 'ar' for n in np.mean(bondTypeList, axis=0)]
            bondData = ['{} {} {} {} \n'.format(i+1,atomid1+1,atomid2+1,bt) for i, ((atomid1,atomid2), bt) in enumerate(zip([(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()], bondType))]
            return bondData

        xyzBlock = self.giveXYZBlock(unit='Angstrom', elementSymbol=True)
        newBondData = generateModifiedBondData(xyzBlock, octetcharge=octetcharge)
        self.__dataDict['BOND'] = newBondData

    def write(self, filePath):
        mol2Data = ''.join([''.join(['@<TRIPOS>{}\n'.format(key)]+value) for key,value in self.__dataDict.items()])
        with open(filePath, mode='w') as f:
            f.write(mol2Data)

