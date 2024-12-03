import re

import numpy as np
from rdkit import Chem

from chemscripts.unit import checkValidUnit, getUnitConversionFactor

class Molecule:
    def __init__(self, atomicnumList=None, symbolList=None, xyzList=None, xyzBlock=None, unit='Angstrom'):
        # Noneチェック
        if unit is None:
            raise ValueError('unit is None')

        table = Chem.GetPeriodicTable()
        if atomicnumList is not None and xyzList is not None:
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
            
        elif symbolList is not None and xyzList is not None:
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
        
        elif xyzBlock is not None:
            if type(xyzBlock) is not str:
                raise TypeError('type of xyzBlock must be str')
            xyzBlock = [line.strip().split() for line in xyzBlock.split()]
            if len(xyzBlock[0]) == 1:
                # ヘッダー行を除去
                xyzBlock = xyzBlock[2:]
            
            if any([len(line)!=4 for line in xyzBlock]):
                raise ValueError('There are rows that does not have 4 columns')
            
            if re.fullmatch('^[0-9]+$', xyzBlock[0][0]):
                # 一列目を原子番号としてパース
                atomicnumList = [int(n) for n,_,_,_ in xyzBlock]
                symbolList = [table.GetElementSymbol(n) for n in atomicnumList]
            else:
                # 一列目を元素記号としてパース
                symbolList = [s for s,_,_,_ in xyzBlock]
                atomicnumList = [table.GetAtomicNumber(s) for s in symbolList]
                
            xyzList = [[float(x),float(y),float(z)] for _,x,y,z in xyzBlock]
            numAtom = len(xyzBlock)

        else:
            raise ValueError('The arguments on a molecular geometry are not specified.')
        
        
        # xyzListチェック
        if type(xyzList) not in [np.ndarray, list, tuple]:
            raise TypeError('type of xyzList must be np.ndarray or list, or tuple')
        xyzArray = np.array(xyzList)
        if xyzArray.dtype.name != 'float64':
            raise TypeError('dtype of xyzList must be float64')
        if len(xyzArray.shape) != 2 or xyzArray.shape[1] != 3:
            raise ValueError('shape of xyzList must be (*,3)')

        # 要素数は一致しているか
        if len(xyzArray) != numAtom:
            raise ValueError('The number of atoms differs between atomicnumList(symbolList) and xyzList')

        if checkValidUnit(unit):
            raise ValueError('Invalid unit: {}'.format(unit))

        # メンバ変数に追加
        self.__numAtom = numAtom
        self.__atomicnumList = atomicnumList
        self.__symbolList = symbolList
        self.__xyzArray = xyzArray
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








