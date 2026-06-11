import re
import copy
import itertools

import numpy as np
from rdkit import Chem

from chemscripts.unit import checkInvalidUnit, getUnitConversionFactor

class Molecule:
    def __init__(self, atomicnumList=None, symbolList=None, xyzList=None, xyzBlock=None, xyzFile=None, charge=None, multiplicity=None, spinState='low', unit='Angstrom'):
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
            if any([type(n) not in [int, np.int32, np.int64] for n in atomicnumList]):
                # 要素は整数のみ
                raise TypeError('type of elements of atomicnumList must be int')
            if any([n<0 for n in atomicnumList]):
                # 要素は非負数のみ
                raise ValueError('elements of atomicnumList must not be negative')
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

        elif xyzBlock is not None or xyzFile is not None:
            if xyzBlock is not None and type(xyzBlock) is not str:
                raise TypeError('type of xyzBlock must be str')
            if xyzFile  is not None and type(xyzFile) is not str:
                raise TypeError('type of xyzFile must be str')
            if xyzFile is None:
                xyzBlock = [line.strip().split() for line in xyzBlock.splitlines()]
            else:
                with open(xyzFile, mode='r') as f:
                    xyzBlock = [line.strip().split() for line in f.readlines()]
            if len(xyzBlock[0]) == 1:
                # 一行目に一つしか値がなかった場合 (= xyzデータではなかった場合)はヘッダー行のあるデータと見做す
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

        if checkInvalidUnit(unit):
            raise ValueError('Invalid unit: {}'.format(unit))

        # charge
        if charge is None:
            charge = 0
        if type(charge) is not int:
            raise TypeError('type of charge must be int')

        # multiplicity, spinState
        if multiplicity is None:
            _spinState = 0
        elif type(multiplicity) is int:
            _spinState = (multiplicity + multiplicity%2)/2 - 1
        else:
            raise TypeError()

        if spinState is None:
            _spinState = 0
        elif type(spinState) is int:
            if spinState < 0:
                raise ValueError()
            _spinState = spinState
        elif type(spinState) is str:
            if re.fullmatch('^low$', spinState):
                _spinState = 0
            elif re.fullmatch('^high[0-9]*$', spinState):
                _spinState = int(re.sub('[^0-9]', '', spinState))
            else:
                raise ValueError()
        else:
            raise ValueError()

        # メンバ変数に追加
        self.__numAtom = numAtom
        self.__atomicnumList = atomicnumList
        self.__symbolList = symbolList
        self.__xyzArray = xyzArray
        self.__charge = charge
        self.__spinState = _spinState
        self.__unit = unit

    def __str__(self):
        return '\n'.join([
                    repr(self),
                    '{} (charge:{},spinState:{})'.format(self.giveStoichiometry(), self.__charge, self.__spinState),
                    ''
                ])

    def giveNumAtom(self):
        return self.__numAtom

    def giveStoichiometry(self):
        """
        give Stoichiometry
        e.g., C3H8O
        """
        from collections import Counter
        count = Counter([s for s,_,_,_ in self.iterateAtoms(unit='Angstrom', elementSymbol=True, atomfilter=None)])

        stoichiometry = ''
        for symb in ['C', 'H', 'N', 'O']:
            if count[symb] == 0:
                continue
            elif count[symb] == 1:
                stoichiometry += '{}'.format(symb)
            else:
                stoichiometry += '{}{}'.format(symb,count[symb])
        for symb, c in count.items():
            if symb not in ['C', 'H', 'N', 'O']:
                stoichiometry += '{}{}'.format(symb,c)
        return stoichiometry

    def iterateAtoms(self, unit='Angstrom', elementSymbol=True, atomfilter=None):
        """
        unit: unit of xyz
        elementSymbol: Whether to output element symbol or atomic number
        atomfilter: Specify element symbols (atomic numbers) to exclude (type is list)
        """
        factor = getUnitConversionFactor(self.__unit, unit)

        xyzlist = (self.__xyzArray * factor).T # shape: (3,n)

        if elementSymbol:
            zipit = zip(self.__symbolList, *xyzlist) # shape: (n,4)
        else:
            zipit = zip(self.__atomicnumList, *xyzlist) # shape: (n,4)

        if atomfilter is None:
            return zipit
        elif type(atomfilter) is list:
            table = Chem.GetPeriodicTable()
            atomfilter = [table.GetAtomicNumber(s) if type(s) is str else s for s in atomfilter] +\
                        [table.GetElementSymbol(n) if type(n) is int else n for n in atomfilter]
            return itertools.filterfalse(lambda v: v[0] in atomfilter, zipit) # shape: (n-m,4)
        else:
            raise TypeError()


    def giveAtomicnumList(self):
        return copy.deepcopy(self.__atomicnumList)

    def giveXYZArray(self, unit='Angstrom', atomfilter=None):
        if atomfilter is None:
            factor = getUnitConversionFactor(self.__unit, unit)
            xyzArray = self.__xyzArray * factor
            return xyzArray

        else:
            return np.array([[x,y,z] for _,x,y,z in self.iterateAtoms(unit=unit, atomfilter=atomfilter)])

    def giveXYZBlock(self, unit='Angstrom', elementSymbol=True, comment='', atomfilter=None, xyzformat=None):
        if xyzformat is None:
            xyzformat = ''
        lineformat = '{} {'+xyzformat+'} {'+xyzformat+'} {'+xyzformat+'}'

        result = [lineformat.format(s,x,y,z) for s, x, y, z in self.iterateAtoms(unit=unit, elementSymbol=elementSymbol, atomfilter=atomfilter)]
        result = [str(len(result)), comment, *result]
        result = '\n'.join(result)

        return result

    def generateRDKitMolObj(self):
        """
        Requires RDKit 2022.09 or higher
        """
        from rdkit.Chem import rdDetermineBonds

        # 不飽和度を計算し、chargeが適切な値か判断
        # 不適切ならば不飽和度に応じて適当に電荷を設定
        #numC = len([n for n in self.__atomicnumList if n in [6]])
        #numH = len([n for n in self.__atomicnumList if n in [1,9,17,35,53]])
        #numN = len([n for n in self.__atomicnumList if n in [7]])
        #numOther = len(self.__atomicnumList) - numC - numH - numN
        #if numOther > 0:
        #    # 不飽和度を計算できない原子が存在するので
        #    # 不飽和度を計算しないで適当に電荷を設定
        #    unsatu = np.nan
        #else:
        #    unsatu =  2 * numC - numH + numN + 2
        #    isHalfint = (unsatu % 2 == 1)
        #    unsatu /= 2

        # Chem.MolFromXYZBlock()は固定小数点しか付けつけないのでformatする
        xyzblock = self.giveXYZBlock(unit='Angstrom', elementSymbol=True, xyzformat=':.8f', atomfilter=[0])
        mol = Chem.MolFromXYZBlock(xyzblock)
        rdDetermineBonds.DetermineBonds(mol, charge=self.__charge)
        return mol

    def generateStandardizedCoordSystem(self, unit='Angstrom', method='PCA'):
        """
        method: PCA, PCA-ignoreHs
        """

        if method == 'PCA':
            from sklearn.decomposition import PCA

            # 原子核の座標を取得
            nucxyz = self.giveXYZArray(unit=unit)
            center = np.mean(nucxyz, axis=0) # shape: (3,)

            # PCA実行
            pca = PCA()
            pca.fit(nucxyz - center)
            # PC3軸目を法線に設定する
            tangent1 = pca.components_[0] # shape: (3,)
            tangent2 = pca.components_[1] # shape: (3,)
            normal = pca.components_[2] # shape: (3,)

        elif method == 'PCA-ignoreHs':
            from sklearn.decomposition import PCA

            # 水素以外の原子核の座標を取得
            nucxyz = self.giveXYZArray(unit=unit, atomfilter=['H'])
            center = np.mean(nucxyz, axis=0) # shape: (3,)

            # PCA実行
            pca = PCA()
            pca.fit(nucxyz - center)
            # PC3軸目を法線に設定する
            tangent1 = pca.components_[0] # shape: (3,)
            tangent2 = pca.components_[1] # shape: (3,)
            normal = pca.components_[2] # shape: (3,)

        else:
            raise ValueError('invalid method')

        return center, tangent1, tangent2, normal



