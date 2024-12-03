import numpy as np
from rdkit import Chem

from chemscripts.unit import checkValidUnit, getUnitConversionFactor

class Molecule:
    def __init__(self, atomicnumList=None, xyzList=None, unit='Angstrom'):
        # Noneチェック
        if atomicnumList is None:
            raise ValueError('atomicnumList is None')
        if xyzList is None:
            raise ValueError('xyzList is None')
        if unit is None:
            raise ValueError('unit is None')

        # atomicnumListチェック
        if type(atomicnumList) in [list, tuple]:
            # listかtupleならnp.ndarrayに変換
            atomicnumList = np.array(atomicnumList)
        if type(atomicnumList) is not np.ndarray:
            raise TypeError('type of atomicnumList must be np.ndarray, list or tuple')
        if atomicnumList.dtype.name not in ['int32','int64'] or any(atomicnumList < 1):
            # 要素は自然数のintのみ
            raise TypeError('elements of atomicnumList must be positive, and the type must be int')
        if len(atomicnumList.shape) != 1:
            raise ValueError('shape of atomicnumList must be (*,)')

        # xyzListチェック
        if type(xyzList) in [list, tuple]:
            # listかtupleならnp.ndarrayに変換
            xyzList = np.array(xyzList)
        if type(xyzList) is not np.ndarray:
            raise TypeError('type of xyzList must be np.ndarray or list, or tuple')
        if xyzList.dtype.name != 'float64':
            raise TypeError('dtype of xyzList must be float64')
        if len(xyzList.shape) != 2 or xyzList.shape[1] != 3:
            raise ValueError('shape of xyzList must be (*,3)')

        # 要素数は一致しているか
        if len(atomicnumList) != len(xyzList):
            raise ValueError('The number of atoms differs between atomicnumList and xyzList')

        if checkValidUnit(unit):
            raise ValueError('Invalid unit: {}'.format(unit))

        # メンバ変数に追加
        self.__numAtom = len(atomicnumList)
        self.__atomicnumArray = atomicnumList
        self.__xyzArray = xyzList
        self.__unit = unit

    def giveNumAtom(self):
        return self.__numAtom

    def iterateAtoms(self, unit='Angstrom', elementSymbol=True):
        factor = getUnitConversionFactor(self.__unit, unit)

        numlist = self.__atomicnumArray # shape: (n,)
        xyzlist = (self.__xyzArray * factor).T # shape: (3,n)

        if elementSymbol:
            table = Chem.GetPeriodicTable()
            # convert atomic number to element symbol
            symblist = [table.GetElementSymbol(int(n)) for n in numlist]

            return zip(symblist, *xyzlist) # shape: (n,4)
        else:
            return zip(numlist, *xyzlist) # shape: (n,4)

    def generateRDKitMolObj(self):
        """
        Requires RDKit 2022.09 or higher
        """
        from rdkit.Chem import rdDetermineBonds

        xyzblock = self.giveXYZBlock(unit='Angstrom', elementSymbol=True)
        mol = Chem.MolFromXYZBlock(xyzblock)
        rdDetermineBonds.DetermineBonds(mol)
        return mol

    def giveXYZBlock(self, unit='Angstrom', elementSymbol=True):
        result = [str(self.__numAtom), 'comment']
        result.extend(
            ['{} {} {} {}'.format(s,x,y,z) for s, x, y, z in self.iterateAtoms(unit=unit, elementSymbol=elementSymbol)]
        )
        result = '\n'.join(result)

        return result








