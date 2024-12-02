import numpy as np

class Molecule:
    def __init__(self, atomicnumList=None, xyzList=None):
        # Noneチェック
        if atomicnumList is None:
            raise ValueError('atomicnumList is None')
        if xyzList is None:
            raise ValueError('xyzList is None')
        
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
        
        # メンバ変数に追加
        self.__numAtom = len(atomicnumList)
        self.__atomicnumList = atomicnumList
        self.__xyzList = xyzList
        
    def giveNumAtom(self):
        return len(self.__atomicnumList)
        
    def iterateAtom(self, unit='Angstrom', elementSymbol=True):
        return None
    
    def generateRDKitMol(self):
        #from rdkit import Chem
        return None
    
    def giveXYZBlock(self, unit='Angstrom'):
        return None
    