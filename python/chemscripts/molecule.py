import numpy as np
from rdkit import Chem

from chemscripts.unit import checkValidUnit, getUnitConversionFactor

class Molecule:
    def __init__(self, atomicnumList=None, symbolList=None, xyzList=None, unit='Angstrom'):
        # Noneチェック
        if atomicnumList is None and symbolList is None:
            raise ValueError('specify either atomicnumList or symbolList')
        if xyzList is None:
            raise ValueError('xyzList is None')
        if unit is None:
            raise ValueError('unit is None')

        table = Chem.GetPeriodicTable()
        if atomicnumList is not None:
            # atomicnumListチェック
            if type(atomicnumList) in [np.ndarray, tuple]:
                # np.ndarrayかtupleならlistに変換
                atomicnumList = list(atomicnumList)
            if type(atomicnumList) is not list:
                raise TypeError('type of atomicnumList must be list, np.ndarray, or tuple')
            if any([type(n) is not int for n in atomicnumList]):
                # 要素は整数のみ
                raise TypeError('type of elements of atomicnumList must be int')
            if any([n<1 for n in atomicnumList]):
                # 要素は自然数のみ
                raise ValueError('elements of atomicnumList must be positive')
            numAtom = len(atomicnumList)
            symbolList = [table.GetElementSymbol(int(n)) for n in atomicnumList]
            
        else:
            # symbolListチェック
            if type(symbolList) in [np.ndarray, tuple]:
                # np.ndarrayかtupleならlistに変換
                symbolList = list(symbolList)
            if type(symbolList) is not list:
                raise TypeError('type of symbolList must be list, np.ndarray, or tuple')
            if any([type(s) is not str for s in symbolList]):
                # 要素は整数のみ
                raise TypeError('type of elements of atomicnumList must be str')
            numAtom = len(symbolList)
            atomicnumList = [table.GetAtomicNumber(s) for s in symbolList]
            
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
        if len(xyzList) != numAtom:
            raise ValueError('The number of atoms differs between atomicnumList(symbolList) and xyzList')

        if checkValidUnit(unit):
            raise ValueError('Invalid unit: {}'.format(unit))

        # メンバ変数に追加
        self.__numAtom = numAtom
        self.__atomicnumList = atomicnumList
        self.__symbolList = symbolList
        self.__xyzArray = xyzList
        self.__unit = unit

    def giveNumAtom(self):
        return self.__numAtom

    def iterateAtoms(self, unit='Angstrom', elementSymbol=True):
        factor = getUnitConversionFactor(self.__unit, unit)

        xyzlist = (self.__xyzArray * factor).T # shape: (3,n)

        if elementSymbol:
            return zip(self.__symbolList, *xyzlist) # shape: (n,4)
        else:
            return zip(self.__atomicnumList, *xyzlist) # shape: (n,4)

    def giveXYZBlock(self, unit='Angstrom', elementSymbol=True, comment=''):
        result = [str(self.__numAtom), comment]
        result.extend(
            ['{} {} {} {}'.format(s,x,y,z) for s, x, y, z in self.iterateAtoms(unit=unit, elementSymbol=elementSymbol)]
        )
        result = '\n'.join(result)

        return result

    def generateRDKitMolObj(self):
        """
        Requires RDKit 2022.09 or higher
        """
        from rdkit.Chem import rdDetermineBonds

        xyzblock = self.giveXYZBlock(unit='Angstrom', elementSymbol=True)
        mol = Chem.MolFromXYZBlock(xyzblock)
        rdDetermineBonds.DetermineBonds(mol)
        return mol








