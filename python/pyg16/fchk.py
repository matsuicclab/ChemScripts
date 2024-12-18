import re
import copy

import numpy as np
from scipy.spatial.distance import squareform

from pyg16.basisfunction import GTOBasis

class Fchk:
    def __init__(self, filePath):
        """
        load fchk file
        """
        # ファイル読み込み
        with open(filePath, mode='r') as f:
            try:
                data = [s.strip('\n') for s in f.readlines()]
            except UnicodeDecodeError as e:
                # binaryではないか
                print('Fchk file, {} may be a binary file'.format(filePath))
                print(e)
                exit()

        # キーワードの出現する行数(index)を取得し、
        # index間毎にdataを分割していく
        keywordLineMatch = [re.match(r'^(.{42}) ([IRC]) (.+)$', line) for line in data]
        keywordLineIndex = [i for i, match in enumerate(keywordLineMatch) if match is not None]
        recordIndexRange = [[start, nextstart] for start,nextstart in zip(keywordLineIndex, keywordLineIndex[1:])]
        recordIndexRange.append([keywordLineIndex[-1],len(data)])

        # [[match object, 分割データ], [*, *], ....]
        recordList = [[keywordLineMatch[start], data[start:nextstart]] for start, nextstart in recordIndexRange]

        # キーワードの値を取得できるように辞書型に変換していく
        def convertRecordToDict(record):
            match = record[0]
            recorddata = record[1]

            keyword = match.group(1).strip(' ')
            valuetype = match.group(2)
            value = re.sub('.+ ', '', match.group(3))

            if len(recorddata) == 1:
                # 単行レコードの場合
                if valuetype == 'R':
                    value = float(value)
                elif valuetype == 'I':
                    value = int(value)

            else:
                # 複数行のレコードの場合
                # valueには要素数(に相当するもの)が入っているため、recorddata[1:]で置換する
                if valuetype == 'C':
                    # レコードタイプがコメント(文字列)だった場合
                    value = ''.join(recorddata[1:])
                else:
                    stringList = ' '.join(recorddata[1:]).split()
                    if valuetype == 'R':
                        value = [float(s) for s in stringList]
                    else: # valuetype == 'I'
                        value = [int(s) for s in stringList]

            return keyword, value

        # [(keyword, dict), (*, *), ....]
        # -> {key1: dict1, key2: dict2, ....}
        recordDict = dict([convertRecordToDict(record) for record in recordList])

        self.__recordDict = recordDict


    def __devideList(self, targetList, ruleList):
        """
        指定された(一次元)リストを複数のリストに分割する
        targetList = [a,b,c,d,e,f,g]
        ruleList = [3,2,1,1]
        のとき、
        [[a,b,c],[d,e],[f],[g]]
        を生成する
        """
        numTotal1 = len(targetList)
        numTotal2 = sum(ruleList)
        if type(numTotal2) not in [int, np.int16, np.int32, np.int64]:
            raise ValueError('ruleList is a list of int')
        if type(ruleList) not in [np.ndarray, list, tuple]:
            raise ValueError('ruleList is a list of int')
        if numTotal1 != numTotal2:
            raise ValueError('The number of elements indicated by targetList and ruleList do not match')

        cumruleList2 = np.cumsum(ruleList)             # array([3, 5, 6, 7] # ruleListが[3,2,1,1]の場合
        cumruleList1 = np.append(0, cumruleList2)[:-1] # array([0, 3, 5, 6])
        result = [targetList[i:j] for i,j in zip(cumruleList1, cumruleList2)]
        return result


    def giveValue(self, key):
        """
        指定されたキーワードに対応する値を返す

        key: キーワード
        return: 数値、または文字列、または数値のリスト
        """
        return self.__recordDict.get(key, None)

    def containsRecord(self, key):
        """
        指定されたキーワードがfchkに含まれているかチェックする

        key: キーワード
        return: boolean
        """
        return key in self.__recordDict.keys()

    def giveRouteSection(self):
        # ルートセクション
        return self.giveValue('Route')

    def giveTitleSection(self):
        # タイトルセクション
        return self.giveValue('Full Title')

    def giveCharge(self):
        # 系全体の電荷
        return self.giveValue('Charge')

    def giveMultiplicity(self):
        # スピン多重度
        return self.giveValue('Multiplicity')

    def giveNumAtoms(self):
        # 全原子数
        return self.giveValue('Number of atoms')

    def giveNumElectrons(self):
        # 全電子数
        # 擬ポテンシャルを張っている場合は内殻電子は含まれないので注意
        return self.giveValue('Number of electrons')

    def giveNumAlphaElectrons(self):
        # 全alpha電子数
        # 擬ポテンシャルを張っている場合は内殻電子は含まれないので注意
        return self.giveValue('Number of alpha electrons')

    def giveNumBetaElectrons(self):
        # 全beta電子数
        # 擬ポテンシャルを張っている場合は内殻電子は含まれないので注意
        return self.giveValue('Number of beta electrons')

    def giveNumBasis(self):
        # 基底関数の総数
        # 擬ポテンシャルの部分はカウントされないので注意
        return self.giveValue('Number of basis functions')

    def isRestrictedOrbital(self):
        # 制限付き計算かどうか
        # beta orbitalが含まれていれば非制限(False)
        return not self.containsRecord('Beta Orbital Energies')

    def giveRotTr(self, toInput=True):
        # input orientationとstandard orientationの間を変換する回転行列と並進ベクトル(bohr単位)
        # return: rot, trans: np.ndarray
        # r' = rot @ r + trans
        value = self.giveValue('RotTr to input orientation')
        if value is None:
            return None

        rot = np.array(value[:9]).reshape(3,3).T
        trans = np.array(value[9:])

        if toInput:
            return rot, trans
        else:
            return rot.T, - rot.T @ trans

    def giveAtomicNums(self):
        # 原子番号リスト
        value = self.giveValue('Atomic numbers')
        return np.array(value)

    def giveNuclearCharges(self):
        # 原子核電荷リスト
        value = self.giveValue('Nuclear charges')
        return np.array(value)

    def giveCoords(self):
        # 原子核座標(bohr単位)
        # Input orientationかStandard orientationか、どちらかを保証することはできない
        value = self.giveValue('Current cartesian coordinates')
        coord = np.array(value).reshape(-1,3)
        return coord

    def giveSCFEnergy(self):
        # SCF energy
        return self.giveValue('SCF Energy')

    def giveTotalEnergy(self):
        # Total energy
        return self.giveValue('Total Energy')

    def giveOrbitalEnergyList(self, merge=False):
        # 各軌道のエネルギーを低い順に返す
        # merge : Trueの場合、alphaとbetaのエネルギーリストを結合
        #       : ['alpha', alpha軌道の何番目の軌道か(int,0始まり), エネルギー]の配列が返る
        #       : Falseの場合、(np.ndarray (alpha, shape: (numBasis,)), np.ndarray (beta, shape: (numBasis,)))が返る
        isRestricted = self.isRestrictedOrbital()

        alphaEnergies = self.giveValue('Alpha Orbital Energies')

        if not isRestricted:
            betaEnergies = self.giveValue('Beta Orbital Energies')
        else:
            betaEnergies = alphaEnergies

        if not merge:
            return alphaEnergies, betaEnergies

        # mergeする場合
        mergedresult = []
        alphai = 0
        betai = 0
        while True:
            alphaEi = alphaEnergies[alphai]
            betaEi = betaEnergies[betai]

            if alphaEi <= betaEi:
                mergedresult.append(['alpha',alphai,alphaEi])
                alphai += 1
            else:
                mergedresult.append(['beta',betai,betaEi])
                betai += 1

            if alphai == len(alphaEnergies):
                mergedresult.append(['beta',betai,betaEi])
                break
            if betai == len(betaEnergies):
                mergedresult.append(['alpha',alphai,alphaEi])
                break

        return mergedresult

    def giveOrbitalCoeffList(self, merge=False):
        """
        各軌道中の基底関数の係数を、軌道エネルギーが低い順に返す
        merge : Trueの場合、alphaとbetaの係数リストを結合
              : ['alpha', alpha軌道の何番目の軌道か(int,0始まり), 軌道係数リスト(np.ndarray)]の配列が返る
              : Falseの場合、(np.ndarray (alpha, shape: (numBasis,numBasis)), np.ndarray (beta, shape: (numBasis,numBasis)))が返る
        """
        isRestricted = self.isRestrictedOrbital()
        numBasis = self.giveNumBasis()

        alphaCoeffs = np.array(self.giveValue('Alpha MO coefficients')).reshape(numBasis, numBasis)

        if not isRestricted:
            betaCoeffs = np.array(self.giveValue('Beta MO coefficients')).reshape(numBasis, numBasis)
        else:
            betaCoeffs = alphaCoeffs

        if not merge:
            return alphaCoeffs, betaCoeffs

        # mergeする場合
        mergedresult = []
        for spin, spinIndex, _ in self.giveOrbitalEnergyList(merge=True):
            if spin == 'alpha':
                coeffs = alphaCoeffs[spinIndex]
            else:
                coeffs = betaCoeffs[spinIndex]
            mergedresult.append([spin, spinIndex, coeffs])

        return mergedresult

    def __mapShellTypeToNumBasis(self,shelltypes):
        def __temp(st):
            if st == 0:
                return 1
            elif st == 1:
                return 3
            elif st == -1:
                return 4
            elif st == 2:
                return 6
            elif st == -2:
                return 5
            elif st == 3:
                return 10
            elif st == -3:
                return 7
            else:
                raise ValueError('shell type:{} is unknown'.format(st))

        return [__temp(x) for x in shelltypes]

    def giveNumBasisEachAtom(self):
        # 各原子にいくつの基底関数が張られているか
        # リストのインデックス0の要素は常に0
        # (インデックスに原子の通し番号を指定できるようにするための処置)

        shelltypes = self.giveValue('Shell types')
        shelltypes = self.__mapShellTypeToNumBasis(shelltypes)

        shellatommap = self.giveValue('Shell to atom map')

        numBasisEachAtomId = [0] * (len(set(shellatommap))+1)
        for atomId, numBasis in zip(shellatommap, shelltypes):
            numBasisEachAtomId[atomId] += numBasis

        return numBasisEachAtomId

    def giveDensityMatrix(self):
        """
        密度行列を返す
        return: np.ndarray (shape: (numBasis, numBasis))
        """
        # densityのリストを取得
        densityList = self.giveValue('Spin SCF Density')
        densityList = np.array(densityList)
        # 三角行列になっているので、元の対称行列に変形
        densityMatrix = np.triu(squareform(densityList))[:-1,1:]
        densityMatrix = densityMatrix + np.tril(densityMatrix.T, k=-1)

        return densityMatrix

    def giveBasisFuncs(self):
        # 基底関数データ取得
        # shell == 同じ指数、同じ核の縮約基底グループ: (1s), (2s 2px 2py 2pz), (3dx2, 3dy2, 3dz2, 3dxy, 3dxz, 3yz), ...
        shelltypes = self.giveValue('Shell types')                         # len == numShell
        numPrimitives = self.giveValue('Number of primitives per shell')   # len == numShell
        exponents = self.giveValue('Primitive exponents')                  # len == numPrimitive
        contractions = self.giveValue('Contraction coefficients')          # len == numPrimitive
        SPcontractions = self.giveValue('P(S=P) Contraction coefficients') # len == numPrimitive
        if SPcontractions is None:
            SPcontractions = [0 for i in contractions]
        coordsshell = np.array(self.giveValue('Coordinates of each shell')).reshape(-1,3) # len == numShell * 3

        # 各shell単位で分割
        exponents = self.__devideList(exponents, numPrimitives)            # len == numShell
        contractions = self.__devideList(contractions, numPrimitives)      # len == numShell
        SPcontractions = self.__devideList(SPcontractions, numPrimitives)  # len == numShell

        basisFuncList = []
        for st, coord, c, spc, ex in zip(shelltypes, coordsshell, contractions, SPcontractions, exponents):
            if st == 0:
                lmnList = [[0,0,0]]
            elif st == 1:
                lmnList = [[1,0,0],[0,1,0],[0,0,1]]
            elif st == -1:
                lmnList = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
            elif st == 2:
                lmnList = [[2,0,0],[0,2,0],[0,0,2],[1,1,0],[1,0,1],[0,1,1]]
            else:
                raise ValueError('not yet implemented')

            if st == -1:
                # contraction修正
                cList = [c, spc, spc, spc]
            else:
                cList = [c for i in range(len(lmnList))]

            for lmn, c in zip(lmnList, cList):
                basisFuncList.append(GTOBasis(coord, lmn, c, ex))

        return basisFuncList

    def calcElectronDensity(self, r):
        """
        指定された座標における電子密度を計算
        r: 電子密度を計算する座標(単位: Bohr): np.ndarray: shape:(*,3) or (3,)
        return: np.ndarray: shape:(*,)
        """
        # 基底関数インスタンスのリストを取得
        basisFuncList = self.giveBasisFuncs()
        # 軌道係数を取得
        numAlphaElec = self.giveNumAlphaElectrons()
        numBetaElec  = self.giveNumBetaElectrons()
        alphaOrbitalCoeffs, betaOrbitalCoeffs = self.giveOrbitalCoeffList(merge=False)
        orbitalCoeffs = np.vstack([alphaOrbitalCoeffs[:numAlphaElec], betaOrbitalCoeffs[:numBetaElec]]) # shape: (numElec, numBasis)

        # 各点で基底関数を評価
        basisFuncValue = np.array([f.calc(r) for f in basisFuncList]) # shape: (numBasis, numPoint)
        # 各点で分子軌道を評価
        moValue = orbitalCoeffs @ basisFuncValue # shape: (numElec, numPoint)
        # 密度計算
        densitydata = np.sum(moValue**2, axis=0) # shape: (numPoint,)

        return densitydata

    def giveElectronDensityCube(self):
        """
        電子密度のcubeデータを生成
        return: Cubeインスタンス
        """
        # 格子点設定
        # TODO
        gridcoord = np.arange(300000).reshape(-1,3)
        densitydata = self.calcElectronDensity(gridcoord)
        cube = Cube()

        return cube

    def giveSpinDensityMatrix(self):
        """
        スピン密度行列を返す
        return: np.ndarray (shape: (numBasis, numBasis))
        """
        # spin densityのリストを取得
        spinDensityList = self.giveValue('Spin SCF Density')
        if spinDensityList is None:
            # 制限付き計算の場合等
            numBasis = self.giveNumBasis()
            spinDensityMatrix = np.zeros([numBasis, numBasis])
            return spinDensityMatrix

        spinDensityList = np.array(spinDensityList)
        # 三角行列になっているので、元の対称行列に変形
        spinDensityMatrix = np.triu(squareform(spinDensityList))[:-1,1:]
        spinDensityMatrix = spinDensityMatrix + np.tril(spinDensityMatrix.T, k=-1)

        return spinDensityMatrix

    def giveMolObject(self, charge=None):
        from rdkit import Chem
        # https://github.com/jensengroup/xyz2mol
        import xyz2mol

        # まず、FCHKからxyz形式を得る
        # 原子番号のリストを取得
        atomicNums = self.giveValue('Atomic numbers')

        # 座標リストを取得
        coords = self.giveValue('Current cartesian coordinates')
        # 単位をa.u.からangstromに変換する
        coords = [0.529177210903 * c for c in coords]
        # [x,y,z]の配列に変換
        coords = [[x,y,z] for x,y,z in zip(coords[0::3],coords[1::3],coords[2::3])]

        # 電荷
        if charge is None:
            charge = self.giveCharge()

        # molオブジェクトを生成
        mols = xyz2mol.xyz2mol(atomicNums, coords, charge=charge)
        if mols is None or len(mols) == 0:
            return None
        elif len(mols) == 1:
            return mols[0]
        else:
            print('multiple mol objects were obtained, return [0]:'+','.join([Chem.MolToSmiles(mol) for mol in mols]))
            return mols[0]


