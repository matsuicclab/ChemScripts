import re
import itertools
import copy

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from chemscripts.molecule import Molecule
from chemscripts.unit import checkInvalidUnit, getUnitConversionFactor

class Cube:
    """
    gaussianのユーティリティcubegenが生成するcubeファイルを読み込むクラス

    """
    def __init__(self, **args):
        if 'filePath' in args.keys():
            self.__init__fromFile(**args)
        elif 'cubeData' in args.keys():
            self.__init__fromCubeData(**args)
        else:
            raise ValueError('args must contain filePath or cubeData')


    def __init__fromFile(self, filePath=None, valueNames=None):
        # 値チェック
        if filePath is None:
            raise ValueError('filePath is None')

        # ファイル読み込み
        with open(filePath, mode='r') as f:
            try:
                data = [s.strip() for s in f.readlines()]
            except UnicodeDecodeError as e:
                # binaryではないか
                print('Cube file, {} may be a binary file'.format(filePath))
                print(e)
                exit()


        # ヘッダー行が存在するか
        # 3行目に'E'が含まれていないかどうかで判別
        if 'E' in data[2]:
            raise IOError('Cube file, {} may not contain header'.format(filePath))



        # 3行目からfloatの配列に変換
        # cubeファイル内の単位はBohrで統一されているので
        # (cubegenでの作成時にBohrとangstromを間接的に指定するがそれは入力パラメータの単位を指している)
        # 読み込み時はunit='Bohr'とすればよく、一部を単位変換することはしない。
        titleData = data[0:2]
        numData = data[2:]
        numData = [s for s in numData if s != '']      # 空行は除去
        numData = [re.split(' +', s) for s in numData] # 数の二重リストに変形する
        numData = [[float(s) for s in l] for l in numData]

        # 結合
        data = titleData + numData

        # ヘッダー行が存在するか(2回目)
        if len(numData[0]) == 6:
            raise IOError('Cube file, {} may not contain header'.format(filePath))

        # 原子数
        numAtom = int(data[2][0])
        # 格子点の基準点（開始位置）
        # [x,y,z]: unit: Bohr
        startingPoint = data[2][1:4]
        # データの次元 (設定されてなければ1にする)
        if len(data[2]) == 4:
            valueDim = 1
        else:
            valueDim = int(data[2][4])

        # 値の名前設定
        # 先に初期化しておく
        self.__valueDim = 0
        self.__valueNames = None
        # チェック
        valueNames =  self.__checkValueNames(valueNames, valueDim)
        self.__valueDim = valueDim
        self.__valueNames = valueNames

        # 格子のステップ数と単位ベクトル
        # [n1,n2,n3]: unit: Bohr
        numGridPoint = [int(l[0]) for l in data[3:6]]
        # [[v1x,v1y,v1z], [v2x,v2y,v2z], [v3x,v3y,v3z]]: unit: Bohr
        stepVector = [l[1:4] for l in data[3:6]]
        # gridインスタンス生成
        self.__cubeGrid = CubeGrid(startingPoint=startingPoint, stepVector=stepVector, numGridPoint=numGridPoint, unit='Bohr')

        # 原子番号、原子座標
        # shape: (numAtom,)
        atomicNumData = [int(l[0]) for l in data[6:6+numAtom]]
        # shape: (numAtom, 3)
        atomXYZData = [[l[2],l[3],l[4]] for l in data[6:6+numAtom]]
        self.__molecule = Molecule(atomicnumList=atomicNumData, xyzList=atomXYZData, unit='Bohr')

        # cubeデータ取り出し
        # 平坦化し、3次元の行列に変換
        self.__cubeData = np.array(
                            list(itertools.chain.from_iterable(data[6+numAtom:]))
                        ).reshape(*self.__numGridPoint,self.__valueDim)


    def __init__fromCubeData(self, cubeGrid=None, cubeData=None, valueDim=1, valueNames=None, moleculeObj=None):
        """
        格子データが既に存在する場合に利用する

        cubeDataはnp.ndarrayであり、
        shapeは1次元(nx*ny*nz*valueDim,)または3次元(nx,ny,nz)(valueDim==1)、
        複数の格子データが存在する場合は4次元(nx,ny,nz,valueDim)とする
        """
        # 値チェック
        if cubeGrid is None:
            raise ValueError('cubeGrid is None')

        if cubeData is None:
            raise ValueError('cubeData is None')

        if moleculeObj is not None and type(moleculeObj) is not Molecule:
            # Noneは許容するが型チェックだけしとく
            raise TypeError('type of moleculeObj must be chemscripts.molecule.Molecule')

        # valueDim
        if type(valueDim) is not int:
            raise TypeError('type of valueDim must be int')
        if valueDim < 1:
            raise ValueError('valueDim must be larger than 0')

        numGridPoint = cubeGrid.giveNumGridPoint()

        # cubeData
        if type(cubeData) not in [np.ndarray, list, tuple]:
            # cubeDataはnp.ndarrayかlistかtupleである必要がある
            raise TypeError('type of cubeData must be np.ndarray or list, or tuple.')
        elif type(cubeData) in [list, tuple]:
            # np.ndarrayに変換
            cubeData = np.array(cubeData)
        if cubeData.dtype.name != 'float64':
            raise TypeError('type of elements of cubeData must be float64')
        if len(cubeData.shape) == 1:
            cubeData = cubeData.reshape(*numGridPoint, valueDim)
        elif len(cubeData.shape) == 3 and valueDim == 1:
            cubeData = cubeData.reshape(*numGridPoint, 1)
        elif len(cubeData.shape) == 4:
            # reshapeする必要はないがnumGridPointに一致するかは確認する
            if cubeData.shape != tuple(list(numGridPoint) + [valueDim]):
                raise ValueError('shape of cubeData must match (nx,ny,nz,valueDim)')
        else:
            # いずれでもなかった場合
            raise ValueError('shape of cubeData must be (nx*ny*nz*valueDim,) or (nx,ny,nz)(valueDim==1), or (nx,ny,nz,valueDim)')

        self.__cubeGrid = cubeGrid
        self.__cubeData = cubeData

        # 値の名前設定
        # 先に初期化しておく
        self.__valueDim = 0
        self.__valueNames = None
        # チェック
        valueNames =  self.__checkValueNames(valueNames, valueDim)
        self.__valueDim = valueDim
        self.__valueNames = valueNames

        self.__molecule = moleculeObj


    def __checkValueNames(self, newValueNames, newValueDim):
        """
        追加するvalueNamesが適切な値になっているかをチェックし、可能ならば修正したもの(intかstrのlist)を返す

        """
        if newValueNames is not None:
            if type(newValueNames) in [int, str]:
                # 型をリストに変換
                newValueNames = [newValueNames]
            if type(newValueNames) in [list, tuple]:
                if len(newValueNames) != newValueDim:
                    raise ValueError('The length of valueName does not fit the number of cube data')
                if any([type(n) not in [int, str] for n in newValueNames]):
                    raise ValueError('ValueNames must be str list or int list.')

            else:
                raise ValueError('ValueNames must be str list or int list.')

        else:
            newValueNames = ['value{}'.format(i) for i in range(self.__valueDim,self.__valueDim+newValueDim)]

        return newValueNames

    def giveCubeData(self):
        """
        cubeデータと値の名前リストのコピーを返す
        """
        return copy.deepcopy(self.__valueNames), copy.deepcopy(self.__cubeData)

    def giveStepVector(self, unit=None):
        """
        cubeデータの実座標復元のための格子ベクトルを返す
        """
        return self.__cubeGrid.giveStepVector(unit=unit)

    def giveNumGridPoint(self):
        """
        グリッドの各方向の点数を返す
        """
        return self.__cubeGrid.giveNumGridPoint()

    def giveStartingPoint(self, unit=None):
        """
        cubeデータの実座標復元のための格子の開始位置を返す
        """
        return self.__cubeGrid.giveStartingPoint(unit=unit)

    def giveNodeCoord(self, unit=None):
        """
        格子座標を返す
        return: np.ndarray (shape: (na, nb, nc, 3))
        """
        return self.__cubeGrid.giveNodeCoord(unit=unit)

    def giveMoleculeObj(self):
        """
        各原子の原子番号と核座標を返す
        return: chemscripts.molecule.Molecule
        """
        return self.__molecule

    def interpolate(self, r, method='linear', unit=None):
        """
        補間値を計算
        https://zenn.dev/tab_ki/articles/interpolation_of_3d_data

        r: 補間する位置ベクトル: np.ndarray (shape: (n, 3))
        return: 補間値: np.ndarray (shape: (n, valueDim))
        """

        # RegularGridInterpolatorを使って補間計算を行う
        # ただし、補間関数を生成するときに、グリッド位置を
        # [0,1,...,na], [0,1,...,nb], [0,1,...,nc]の3配列で指定する仕様であるため、
        # 実格子座標(r)を直交系座標(p)へ変換する(__convertCoordToOrtho)

        # 補間関数の生成
        na, nb, nc = self.giveNumGridPoint()
        # bounds_error: Falseの場合、補外することになった場合適当な値を代わりにセットする(デフォルト: np.nan)
        interp = RegularGridInterpolator((np.arange(na), np.arange(nb), np.arange(nc)), self.__cubeData, bounds_error=False)

        # 指定された座標を直交系に変換
        p = self.__cubeGrid.convertRealCoordToGridCoord(r, unit=unit)

        # 補間手法を設定
        interp.method = method
        # 補間値を計算
        return interp(p)


    def addCubeData(self, newCubeData, newValueNames=None):
        """
        cubeデータを新しく追加する
        newCubeData: numpy配列((x方向の点数)*(y方向の点数)*(z方向の点数)*(値の数)の次元)
        """
        # numpy配列か
        if type(newCubeData) is not np.ndarray:
            raise ValueError('Type of newCubeData must be numpy array')

        # 次元をチェック
        if len(self.__cubeData.shape) != len(newCubeData.shape) or self.__cubeData.shape[:-1] != newCubeData.shape[:-1]:
            raise ValueError('Does not fit the xyz dimension of the existing cubeData.')

        # 値の数を取得
        newValueDim = newCubeData.shape[-1]

        # 値の名前設定
        # チェック
        newValueNames = self.__checkValueNames(newValueNames, newValueDim)
        self.__valueDim += newValueDim
        self.__valueNames.extend(newValueNames)

        # 問題なければcubeDataに追加
        self.__cubeData = np.block([self.__cubeData, newCubeData])

    def write(self, cubeFilePath, header=True):
        """
        CubeデータをCubeファイルに書き出し

        cubeFilePath: 書き出し先
        header: ヘッダーを入れるかどうか
        """
        if type(header) is not bool:
            raise ValueError('type of \'header\' must be bool')
        if type(cubeFilePath) is not str:
            raise ValueError('type of \'cubeFilePath\' must be str')

        """
        Cubeファイル仕様(https://qiita.com/sahayan/items/ac700c5478920e2cd3dd)
        - 1,2行目...タイトル(コメント)
        - 3行目...{原子数} {格子点の(0,0,0)位置x} {y} {z} {Cubeデータの次元数}
        - 4行目...{1つ目の格子ベクトル方向への格子点数} {1つ目の格子ベクトルx} {y} {z}
        - 5行目...2つ目の格子ベクトル
        - 6行目...3つ目の格子ベクトル
        - 7行目~...{原子番号} {価電子数} {x} {y} {z}
            - 原子データが存在しない場合は書き出さないようにする
            - ここまでヘッダー
        - 以降...{Cubeデータ(0,0,0,0)} {(0,0,0,1)} {(0,0,0,2)} {(0,0,0,3)} {(0,0,0,4)} {(0,0,0,5)}
            - 各行6個ずつ値を出力
            - (0,0,n3,d)まで出力したら改行して(0,1,0,0)から再開
        """
        if header:
            header1 = 'Cube Data generated by chemscripts.pyg16.cube.Cube.write()\n'
            header2 = 'value:{}\n'.format(self.__valueNames)
            header3 = '{} {} {} {} {}\n'.format(self.__molecule.giveNumAtom(), *self.giveStartingPoint(unit='Bohr'), self.__valueDim)
            header456 = ''.join(['{} {} {} {}\n'.format(n,*v) for n,v in zip(self.giveNumGridPoint(),self.giveStepVector(unit='Bohr'))])
            if self.__atomicNumData is not None:
                header7 = ''.join(['{} {} {} {} {}\n'.format(n,n,x,y,z) for n,x,y,z in self.__molecule.iterateAtoms(unit='Bohr',elementSymbol=False)])
            else:
                header7 = ''

            headerContents = ''.join([header1,header2,header3,header456,header7])
        else:
            headerContents = ''

        # format用のテンプレートを予め最後まで作っておいて、一気にcubeDataを流し込む感じで作る
        # (0,0,n3,d)までの間に6個値を出力する行の数と、端数の行の数(0 or 1)を取得し、(0,0,n3,d)までのテンプレートをまず作成
        numValueTon3d = self.__cubeData.shape[2] * self.__cubeData.shape[3]
        numLineSixValue = numValueTon3d // 6
        numFraction = numValueTon3d % 6 # 端数行に含まれる値の数
        numLineFraction = np.sign(numFraction) # 端数行の数(0 or 1)
        sixValueLineTemplate = ''.join(['{} ' for _ in range(6)]) + '\n'
        fractionLineTemplate = ''.join(['{} ' for _ in range(numFraction)]) + '\n'
        formatTemplate = ''.join([sixValueLineTemplate for _ in range(numLineSixValue)] + [fractionLineTemplate for _ in range(numLineFraction)])
        formatTemplate = ''.join(itertools.repeat(formatTemplate, self.__cubeData.shape[0] * self.__cubeData.shape[1]))

        mainContents = formatTemplate.format(*self.__cubeData.reshape(-1))

        # 書き出し
        with open(cubeFilePath, mode='w') as f:
            f.write(headerContents)
            f.write(mainContents)



class CubeGrid:
    def __init__(self, **args):
        if 'moleculeObj' in args.keys():
            self.__init__fromMoleculeObj(**args)
        elif '--' in args.keys():
            self.__init__fromParam(**args)
        else:
            raise ValueError('args must contain moleculeObj or --')
    
    def __init__fromMoleculeObj(self, moleculeObj=None, axesMethod=None, step=0.5, padding=3, unit='Angstrom'):
        # molecule
        if moleculeObj is None:
            raise ValueError()
        if type(moleculeObj) is not Molecule:
            raise TypeError()
        # unit
        if checkInvalidUnit(unit):
            raise ValueError('Invalid unit: {}'.format(unit))
        
        # stepVector決定
        if axesMethod is None or axesMethod == 'Direct':
            v1,v2,v3 = np.diag([1,1,1])
            
        elif axesMethod == 'PCA' or axesMethod == 'PCA-ignoreHs':
            _, v1, v2, v3 = moleculeObj.generateStandardizedCoordSystem(method=axesMethod)
            
        else:
            raise ValueError()
        stepVector = np.array([v1,v2,v3])

        # startingPoint決定
        # 各原子のxyz座標を取得
        xyzAtom = moleculeObj.giveXYZArray(unit=unit)
        # グリッド座標(原点位置はxyz座標と共通)を取得
        V = stepVector.T
        Vinv = np.linalg.inv(V)
        ijkAtom = (Vinv @ xyzAtom.T).T # shape: (n, 3)
        # グリッド座標で最も端の座標を取得
        minijkAtom = np.min(ijkAtom, axis=0) # shape: (3,)
        maxijkAtom = np.max(ijkAtom, axis=0) # shape: (3,)
        # 端の位置のxyz座標を取得
        startingPoint = (V @ minijkAtom.T).T - padding * (v1 + v2 + v3)
        endingPoint = (V @ maxijkAtom.T).T + padding * (v1 + v2 + v3)
        
        self.__init__fromParam(startingPoint=startingPoint, stepVector=stepVector, endingPoint=endingPoint, unit=unit)
        
        # TODO stepVectorの単位に関する解釈がおかしい気がする
        # 恐らくSliceクラスでのnormalVectorなどでは無次元なのに対して、Cubeでは有次元なのが混乱の原因かと
        # 引数stepが未使用

    
    def __init__fromParam(self, startingPoint=None, stepVector=None, numGridPoint=None, endingPoint=None, unit=None):

        # startingPoint
        # None判定
        if startingPoint is None:
            raise ValueError('startingPoint is None')
        # 型チェック
        if type(startingPoint) not in [np.ndarray, list, tuple]:
            # startingPointはnp.ndarrayかlistかtupleである必要がある
            raise TypeError('type of startingPoint must be np.ndarray or list, or tuple.')
        elif type(startingPoint) in [list, tuple]:
            # np.ndarrayに変換
            startingPoint = np.array(startingPoint)
        if startingPoint.shape != (3,):
            raise ValueError('shape of startingPoint must be (3,)')
        if startingPoint.dtype.name not in ['float64', 'int32', 'int64']:
            raise TypeError('dtype of startingPoint must be float64 or int32')


        nonetuple = (stepVector is None, numGridPoint is None, endingPoint is None)
        if nonetuple == (True, True, True) or nonetuple == (True, True, False) or nonetuple == (True, False, True) or nonetuple == (False, True, True):
            raise ValueError('At least two of stepVector, numGridPoint, and endingPoint must be specified')

        # stepVector
        if stepVector is not None:
            # 型チェック
            if type(stepVector) not in [np.ndarray, list, tuple]:
                # stepVectorはnp.ndarrayかlistかtupleである必要がある
                raise TypeError('type of stepVector must be np.ndarray or list, or tuple.')
            elif type(stepVector) in [list, tuple]:
                # np.ndarrayに変換
                stepVector = np.array(stepVector)
            if stepVector.shape != (3,3):
                raise ValueError('shape of stepVector must be (3,3)')
            if stepVector.dtype.name not in ['float64', 'int32', 'int64']:
                raise TypeError('dtype of stepVector must be float64 or int32, or int64')


        # numGridPoint
        if numGridPoint is not None:
            # 型チェック
            if type(numGridPoint) is not np.ndarray:
                # listかtupleならnp.ndarrayに変換
                if type(numGridPoint) not in [list, tuple]:
                    raise TypeError('type of numGridPoint must be np.ndarray, list or tuple')
                else:
                    numGridPoint = np.array(numGridPoint)
            if numGridPoint.shape != (3,):
                raise ValueError('length of numGridPoint must be 3')
            if numGridPoint.dtype.name not in ['int32', 'int64']:
                raise TypeError('dtype of numGridPoint must be int32 or int64')
            if any(numGridPoint < 1):
                raise ValueError('the number of grids in one direction must be greater than or equal to 1')

        # endingPoint
        if endingPoint is not None:
            # 型チェック
            if type(endingPoint) not in [np.ndarray, list, tuple]:
                # endingPointはnp.ndarrayかlistかtupleである必要がある
                raise TypeError('type of endingPoint must be np.ndarray or list, or tuple.')
            elif type(endingPoint) in [list, tuple]:
                # np.ndarrayに変換
                endingPoint = np.array(endingPoint)
            if endingPoint.shape != (3,):
                raise ValueError('shape of endingPoint must be (3,)')
            if endingPoint.dtype.name not in ['float64', 'int32', 'int64']:
                raise TypeError('dtype of endingPoint must be float64 or int32')

        # 欠けている量を決定
        if nonetuple == (True, False, False):
            # stepVector決定
            if any((endingPoint-startingPoint)[np.where(numGridPoint==1)]==0):
                raise ValueError('')
            stepVector = np.diag(endingPoint-startingPoint)/(numGridPoint-1) # shape: (3,3): [v1,v2,v3]

        elif nonetuple == (False, True, False):
            # numGridPoint決定
            numGridPoint = np.linalg.inv(stepVector.T)@(endingPoint-startingPoint) + 1
            # 切り上げ
            numGridPoint = np.ceil(numGridPoint).astype('int64')
            # endingPoint更新
            endingPoint = startingPoint + stepVector.T @ (numGridPoint-1) # r1 = r0 + [v1,v2,v3]@[n1,n2,n3]

        else:
            # if nonetuple == (False, False, True) or nonetuple == (False, False, False):
            # endingPointは最終的に計算結果に影響しないが一応決定する
            # endingPointも含めて全て指定されていた場合はendingPointを上書きする
            endingPoint = startingPoint + stepVector.T @ (numGridPoint-1)

        # unit
        if checkInvalidUnit(unit):
            raise ValueError('Invalid unit: {}'.format(unit))

        # メンバ変数に追加
        self.__startingPoint = startingPoint
        self.__endingPoint = endingPoint
        self.__stepVector = stepVector
        self.__numGridPoint = numGridPoint
        self.__unit = unit


    def giveStartingPoint(self, unit=None):
        """
        cubeデータの実座標復元のための格子の開始位置を返す
        """
        factor = getUnitConversionFactor(oldunit=self.__unit, newunit=unit)
        return copy.deepcopy(self.__startingPoint*factor)

    def giveEndingPoint(self, unit=None):
        """
        cubeデータの実座標復元のための格子の終了位置を返す
        """
        factor = getUnitConversionFactor(oldunit=self.__unit, newunit=unit)
        return copy.deepcopy(self.__endingPoint*factor)

    def giveStepVector(self, unit=None):
        """
        cubeデータの実座標復元のための格子ベクトルを返す
        """
        factor = getUnitConversionFactor(oldunit=self.__unit, newunit=unit)
        return copy.deepcopy(self.__stepVector*factor)

    def giveNumGridPoint(self):
        """
        グリッドの各方向の点数を返す
        """
        return copy.deepcopy(self.__numGridPoint)

    def giveNodeCoord(self,unit=None):
        """
        格子座標を返す
        return: np.ndarray (shape: (na, nb, nc, 3))
        """
        na, nb, nc = self.__numGridPoint
        a,b,c = np.meshgrid(np.arange(na),np.arange(nb),np.arange(nc), indexing='ij')
        sv1 = self.__stepVector[0]
        sv2 = self.__stepVector[1]
        sv3 = self.__stepVector[2]
        nodeCoord = self.__startingPoint + a[:,:,:,np.newaxis] * sv1 + b[:,:,:,np.newaxis] * sv2 + c[:,:,:,np.newaxis] * sv3

        factor = getUnitConversionFactor(oldunit=self.__unit, newunit=unit)
        return nodeCoord * factor

    def convertRealCoordToGridCoord(self, r, unit=None):
        """
        実座標(xyz)をグリッド座標(ijk)へ変換する

        r: np.ndarray (shape: (n, 3) or (3,))
        return: np.ndarray (shape: (n, 3) or (3,))
        """
        # -> r = r0 + [v1,v2,v3] @ p
        #    p = V^-1 (r-r0)
        factor = getUnitConversionFactor(oldunit=unit, newunit=self.__unit)
        r = r * factor

        if r.shape == (3,):
            _r = r.reshape(1,3)
        else:
            _r = r

        V = self.__stepVector.T # stepVectorは[v1,v2,v3]になっているが、v1等は横ベクトルになっているので縦ベクトルに転置
        Vinv = np.linalg.inv(V)
        p = (Vinv @ (_r - self.__startingPoint).T).T # shape: (n, 3)

        if r.shape == (3,):
            return p.reshape(3)
        else:
            return p

    def convertGridCoordToRealCoord(self, p, unit=None):
        """
        グリッド座標(ijk)を実座標(xyz)へ変換する (逆変換)

        p: np.ndarray (shape: (n, 3) or (3,))
        return: np.ndarray (shape: (n, 3) or (3,))
        """
        factor = getUnitConversionFactor(oldunit=self.__unit, newunit=unit)

        if p.shape == (3,):
            _p = p.reshape(1,3)
        else:
            _p = p

        V = self.__stepVector.T
        r = self.__startingPoint + (V@_p.T).T # shape: (n, 3)

        r = r * factor

        if p.shape == (3,):
            return r.reshape(3)
        else:
            return r


class Slice:
    """
    Cubeから生成されたスライス面のデータクラス
    """
    def __init__(self, cube, pos=None, normal=None, pcaAuto=False, unit=None):
        """
        スライスデータを生成
        pos: np.ndarray or list (shape: (3,))
        normal: np.ndarray or list (shape: (3,))
        pcaAuto: 原子座標に関してPCA分析を使ってスライス面を自動決定

        """
        if type(pos) is list:
            pos = np.array(pos)
        if type(normal) is list:
            normal = np.array(normal)

        # posとnormalを設定
        if pcaAuto:
            if unit is None:
                unit = 'Angstrom'

            # 原子核の座標を取得
            molecule = cube.giveMoleculeObj()
            pos, tanVector, tan2Vector, normVector = molecule.generateStandardizedCoordSystem(unit=unit,method='PCA-ignoreHs')
            
        else:
            if pos is None or normal is None:
                raise ValueError()
            elif type(pos) is not np.ndarray or type(normal) is not np.ndarray:
                raise TypeError()
            elif pos.shape != (3,) or normal.shape != (3,):
                raise ValueError()

            # 面の接線ベクトルを生成
            normVector = normal / np.linalg.norm(normal)
            for tanVector in np.diag([1,1,1]): # 接線ベクトルの候補3つを順に判定
                # 直交化
                tanVector = tanVector - (tanVector@normVector) * normVector
                if np.all(tanVector == 0):
                    # もしもtanVectorとnormVectorが線形従属ならば
                    continue
                else:
                    # 線形独立の場合
                    break
            tan2Vector = np.cross(normVector, tanVector)
    
            # 規格化
            tanVector  = tanVector  / np.linalg.norm(tanVector)
            tan2Vector = tan2Vector / np.linalg.norm(tan2Vector)
            
        self.__tanVector = tanVector
        self.__tan2Vector = tan2Vector
        self.__normVector = normVector

        # スライスの座標を生成
        #
        # スライスの幅がcubeの範囲を上回るように調整するために、
        # 1. cubeの格子辺12本の端点v1,v2を用意
        # 2. 格子辺12本とスライスの交点を計算
        #     pos + alpha * t1 + beta * t2 = v1 + gamma * (v2 - v1)
        #     -> [t1, t2, v1-v2]@[alpha,beta,gamma] = v1 - pos
        # 3. 内分点(0<=gamma<=1)になるものだけを抽出
        # 4. 内分点の座標からスライスの幅を決める
        # 5. スライスの幅から適当な間隔でスライス上の座標を生成
        # 6. 生成した座標がcubeの範囲内に存在するか確認

        # 1.
        #
        # nsv1 = a * sv1 (== v1-v2)
        # -> nsv2, nsv3を使ってv1を生成([0,nsv2,nsv3,nsv2+nsv3])、そこからnsv1を足してv2を生成
        svs = cube.giveStepVector(unit=unit) # shape: (3,3)
        nsvs = cube.giveNumGridPoint()[:,np.newaxis] * svs # shape: (3,3)
        arr = nsvs[[1,2, 0,2, 0,1]].reshape(3,2,3) # shape: (3,2,3)
        arr2 = np.sum(arr, keepdims=True, axis=1)  # shape: (3,1,3): nsv2+nsv3に相当
        v1s = cube.giveStartingPoint(unit=unit) + np.hstack([np.zeros(arr2.shape), arr, arr2]) # shape: (3,4,3)
        v2s = nsvs[[0, 1, 2],np.newaxis,:] + v1s           # shape: (3,4,3)

        # 2.
        #
        # 交点座標計算
        # 12本の交点を同時に計算
        delta = (v2s-v1s).reshape(-1,3) # shape: (12,3)
        A = np.stack([np.tile(tanVector, (12,1)), np.tile(tan2Vector, (12,1)), -delta], axis=2) # shape: (12,3,3)
        b = (v1s.reshape(-1,3) - pos)[:,:,np.newaxis] # shape: (12,3,1)
        Apinv = np.linalg.pinv(A) # shape: (12,3,3)
        abg = Apinv @ b # alpha,beta,gamma, shape: (12,3,1)
        # gammaの情報が重要なので、そこだけ取り出す
        gamma = abg[:,2,:].reshape(-1)

        # 3.
        #
        # 内分点に絞り込み
        # Aが正則でない場合はgammaがinfになるように修正する
        gamma = np.where(np.abs(np.sign(np.linalg.det(A)))==0,np.inf,gamma)
        isInternal = (gamma>=0) & (gamma<=1)
        # 交点を実際に計算する
        gamma = gamma[isInternal] # shape: (*,3)
        v1s = v1s.reshape(-1,3)[isInternal] # shape: (*,3)
        delta = delta[isInternal] # shape: (*,3)
        intersec = v1s + gamma.reshape(-1,1) * delta # shape: (*,3)

        # 4.
        #
        # 交点の中心centerから各交点を見て、t1,t2方向の幅を決める。
        center = np.mean(intersec, axis=0) # shape: (3,)
        self.__center = center
        #
        intersec_c1, intersec_c2 = self.project3DCoordToSlice(intersec) # shape: (*,)
        intersec_c = np.stack([intersec_c1, intersec_c2]).T # shape: (*,2)
        mincomponent = np.min(intersec_c, axis=0) # shape: (2,)
        maxcomponent = np.max(intersec_c, axis=0) # shape: (2,)

        # 5.
        #
        # min~maxの間を適当な間隔で生成
        numpatch = 100
        c1, c2 = np.meshgrid(np.linspace(mincomponent[0], maxcomponent[0], numpatch), np.linspace(mincomponent[1], maxcomponent[1], numpatch)) # shape: (numpatch, numpatch)
        self.__numpatch = numpatch
        self.__c1 = c1
        self.__c2 = c2
        # スライス上の点の座標を生成
        r = self.convert2DCoordTo3DCoord(c1, c2).reshape(numpatch,numpatch,3) # shape: (numpatch,numpatch,3)

        # 6.
        #
        # スライス上の値を計算し、Sliceを生成
        value = cube.interpolate(r.reshape(-1,3),unit=unit).reshape(numpatch, numpatch)

        self.__r = r
        self.__value = value
        self.__unit = unit

    def give3DCoord(self,unit=None):
        """
        return: shape: (numpatch,numpatch,3)
        """
        factor = getUnitConversionFactor(self.__unit, unit)
        return copy.deepcopy(self.__r * factor)
    def give2DCoord(self,unit=None):
        """
        return: (shape: (numpatch,numpatch), shape: (numpatch,numpatch))
        """
        factor = getUnitConversionFactor(self.__unit, unit)
        return copy.deepcopy(self.__c1*factor), copy.deepcopy(self.__c2*factor)
    def convert2DCoordTo3DCoord(self, c1, c2):
        """
        c1, c2: shape: (*,)
        return: shape: (*,3)
        """
        r = self.__center + c1.reshape(-1,1) * self.__tanVector + c2.reshape(-1,1) * self.__tan2Vector
        return r
    def project3DCoordToSlice(self, r):
        """
        r: shape: (*,3)
        return: c1,c2 : shape: (*,)
        """
        components = (r - self.__center) @ np.stack([self.__tanVector, self.__tan2Vector], axis=1) # shape: (*,2)
        return components[:,0], components[:,1]

    def giveSliceValue(self):
        """
        return: shape: (numpatch,numpatch)
        """
        return copy.deepcopy(self.__value)
    def give3DRange(self,unit=None):
        """
        return: shape: (3,2) : [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
        """
        factor = getUnitConversionFactor(self.__unit, unit)

        _value = self.__value.reshape(-1) # shape: (numpatch**2,)
        _r = self.__r.reshape(-1,3)       # shape: (numpatch**2,3)
        _r_notnan = _r[~np.isnan(_value)]
        result = np.vstack([np.min(_r, axis=0), np.max(_r, axis=0)]) # shape: (2,3)
        return result.T * factor # shape: (3,2)
    def giveSliceCenter(self,unit=None):
        """
        return: shape: (3,)
        """
        factor = getUnitConversionFactor(self.__unit, unit)
        return copy.deepcopy(self.__center * factor)
    def giveSliceNormalVector(self):
        """
        normal vector (normalized)
        return: shape: (3,)
        """
        return copy.deepcopy(self.__normVector)
    def giveSliceTangentVector(self):
        """
        tangent vector (normalized)
        return: shape: (3,)
        """
        return copy.deepcopy(self.__tanVector)


class CubeVisualizer:
    """
    仮称
    """
    def __init__(self):
        pass

    def __giveRadius(self, n, scale,unit=None):
        """
        原子半径
        周期      2r[A]
        1(1-2)       0.7878
        2(3-10)      1.0918
        3(11-18)     1.1801
        4(19-36)     1.2477
        5(37-54)     1.29
        6(55-86)     1.34
        7(87-118)    1.3903


        n: 原子番号 (np.ndarray, shape:(n,))
        """
        def convertToRadius(n):
            if n<=2:
                return 0.7878/2
            elif n<=10:
                return 1.0918/2
            elif n<=18:
                return 1.1801/2
            elif n<=36:
                return 1.2477/2
            elif n<=54:
                return 1.29/2
            elif n<=86:
                return 1.34/2
            else:
                return 1.3903/2
        factor = getUnitConversionFactor('Angstrom', unit)

        return np.array([convertToRadius(ni) for ni in n]) * scale * factor

    def __giveColorConfig(self, n):
        """
        原子番号から色を設定
        """
        colorDict = {
            1 : [204,204,204],
            2 : [217,255,255],
            #
            3 : [204,125,255],
            4 : [204,255,  0],
            5 : [255,181,181],
            6 : [143,143,143],
            7 : [ 25, 25,229],
            8 : [229,  0,  0],
            9 : [178,255,255],
            10: [176,227,245],
            #
            11: [171, 92,242],
            12: [178,204,  0],
            13: [209,166,166],
            14: [127,153,153],
            15: [255,127,  0],
            16: [255,199, 41],
            17: [ 25,240, 25],
            18: [127,209,227],
            #
            19: [143, 64,212],
            20: [153,153,  0],
            21: [229,229,227],
            22: [191,194,199],
            23: [166,166,171],
            24: [138,153,199],
            25: [156,122,199],
            26: [127,122,199],
            27: [ 92,110,255],
            28: [ 92,122,194],
            29: [255,122, 97],
            30: [125,127,176],
            31: [194,143,143],
            32: [102,143,143],
            33: [189,127,227],
            34: [255,161,  0],
            35: [166, 33, 33],
            36: [ 92,186,209],
            #
            37: [112, 46,176],
            38: [127,102,  0],
            39: [148,252,255],
            40: [148,224,224],
            41: [115,194,201],
            42: [ 84,181,181],
            43: [ 59,158,168],
            44: [ 36,143,150],
            45: [ 10,125,140],
            46: [  0,105,133],
            47: [153,199,255],
            48: [255,217,143],
            49: [166,117,115],
            50: [102,127,127],
            51: [158, 99,181],
            52: [212,122,  0],
            53: [148,  0,148],
            54: [ 66,158,176],
            #
            55: [ 87, 23,143],
            56: [102, 51,  0],
            57: [112,222,255],
            58: [255,255,199],
            59: [217,255,199],
            60: [199,255,199],
            61: [163,255,199],
            62: [143,255,199],
            63: [ 97,255,199],
            64: [ 69,255,199],
            65: [ 48,255,199],
            66: [ 31,255,181],
            67: [  0,255,181],
            68: [  0,229,117],
            69: [  0,212, 82],
            70: [  0,191, 56],
            71: [  0,171, 36],
            72: [ 76,194,255],
            73: [ 76,166,255],
            74: [ 38,148,214],
            75: [ 38,125,171],
            76: [ 38,102,150],
            77: [ 23, 84,135],
            78: [ 23, 92,143],
            79: [255,209, 36],
            80: [181,181,194],
            81: [166, 84, 76],
            82: [ 87, 89, 97],
            83: [158, 79,181],
            84: [171, 92,  0],
            85: [117, 79, 69],
            86: [ 66,130,150],
            #
            87: [ 66,  0,102],
            88: [ 76, 25,  0],
            89: [112,171,250],
            90: [  0,186,255],
            91: [  0,161,255],
            92: [  0,143,255],
            93: [  0,127,242],
            94: [  0,107,242],
            95: [ 84, 92,242],
            96: [120, 92,227],
            97: [138, 94,227],
            98: [161, 54,212],
            99: [168, 43,199],
            100:[178, 31,186],
            101:[178, 13,166],
            102:[189, 13,135],
            103:[199,  0,102],
            104:[255,127,127],
            105:[229,102,102],
            106:[204, 76, 76],
            107:[178, 51, 51],
            108:[153, 25, 25],
            109:[140,  0,  0],
            110:[127,  0,  0],
            111:[115,  0,  0]#,
            #
            #112: [,,],
            #113: [,,],
            #114: [,,],
            #115: [,,],
            #116: [,,],
            #117: [,,],
            #118: [,,],
            #119: [,,],
        }
        colorConfig = ['rgb({},{},{})'.format(*colorDict[ni]) for ni in n]

        return colorConfig

    def giveMoleculeSurfacePlot(self, cube, scale=0.75, unit=None):
        """
        分子の描画オブジェクトを生成
        scale: 原子半径のスケールを調整

        return: list of plotly.graph_objects.Surface
        """
        import plotly.graph_objects as go

        molecule = cube.giveMoleculeObj()
        atomicnumList = molecule.giveAtomicnumList()
        centerList = molecule.giveXYZArray(unit=unit)
        radiusList = self.__giveRadius(atomicnumList, scale)
        colorList = self.__giveColorConfig(atomicnumList)

        datList = []
        for c, r, rgb in zip(centerList, radiusList, colorList):
            theta = np.linspace(0,np.pi,15)
            phi = np.linspace(0,2*np.pi,15)
            theta, phi = np.meshgrid(theta, phi)
            x = r * np.sin(theta) * np.cos(phi) + c[0]
            y = r * np.sin(theta) * np.sin(phi) + c[1]
            z = r * np.cos(theta) + c[2]

            surfaceDat = go.Surface(
                x=x, y=y, z=z,
                colorscale=[[0,rgb],[1,rgb]],
                showscale=False
            )
            datList.append(surfaceDat)

        return datList

    def giveIsosurfacePlot(self, cube, value, enableNegative=True, unit=None):
        """
        等値面プロットを生成
        return: plotly.graph_objects.Isosurface
        """
        import plotly.graph_objects as go

        # 各節点の座標を取得
        node = cube.giveNodeCoord(unit=unit)
        # 各点における値を取得
        _, cubedata = cube.giveCubeData()

        if enableNegative:
            isomax = abs(value)
            isomin = -abs(value)
        else:
            isomax = isomin = value

        cubeIsosurfaceDat = go.Isosurface(
            x=node[:,:,:,0].reshape(-1),
            y=node[:,:,:,1].reshape(-1),
            z=node[:,:,:,2].reshape(-1),
            value=cubedata[:,:,:,0].reshape(-1),
            isomin=isomin,
            isomax=isomax,
            caps=dict(x_show=False, y_show=False, z_show=False)
        )

        return cubeIsosurfaceDat

    def giveSlicePlot(self, slice, unit=None):
        """
        スライスプロットを生成

        return: plotly.graph_objects.Surface
        """
        import plotly.graph_objects as go

        r = slice.give3DCoord(unit=unit)
        value = slice.giveSliceValue()

        sliceDat = go.Surface(
            x=r[:,:,0], y=r[:,:,1], z=r[:,:,2],
            surfacecolor=value
        )

        return sliceDat

    def giveIsolinesPlot(self, slice, numIsoline=None, stepIsoline=1, cutIsolineNote=None, thresholdNoteArrow=-np.inf, unit=None):
        """
        スライス上の等値線プロットを生成
        numIsoline: 等値線の数を設定
        stepIsoline: 等値線の間隔を設定
        cutIsolineNote: 等値線のレベル表記を省く対象を関数で設定
        thresholdNoteArrow: 等値線のレベル表記に矢印を用いる閾値を設定(閾値以下に対して矢印を用いる)

        return: list of plotly.graph_objects.Scatter3d, list of annotations
        """
        import plotly.graph_objects as go
        import matplotlib.pyplot as plt

        if cutIsolineNote is None:
            cutIsolineNote = lambda level : False
        if thresholdNoteArrow is None:
            thresholdNoteArrow = -np.inf

        # スライスのデータを取得
        c1, c2 = slice.give2DCoord(unit=unit)
        value = slice.giveSliceValue()
        # スライスの代表長さ
        sliceRepLength1 = np.max(c1)-np.min(c1)
        sliceRepLength2 = np.max(c2)-np.min(c2)

        # isolineのレベルを決定
        level_min = np.nanmin(value)
        level_max = np.nanmax(value)
        if numIsoline is not None:
            # numIsolineが指定されている場合はレベルの数を指定する
            levels = np.linspace(level_min, level_max, numIsoline)
        else:
            # stepIsolineが指定されている場合はレベルの間隔を指定する
            # levelがstepIsolineの整数倍になるように調整
            levels = np.arange(
                        np.ceil(level_min/stepIsoline)*stepIsoline,
                        (np.floor(level_max/stepIsoline)+1)*stepIsoline,
                        stepIsoline
                    )

        #
        # matplotのcontourを使って座標を取得する
        isolines = plt.contour(c1, c2, value, levels=levels)
        plt.close() # これがないとplt.contour()が描画されるのでそれを防ぐ
        # 等値線の座標を折れ線グラフで描画する
        isolineDatList = []
        annotationList = []
        for i in range(len(isolines.collections)):
            paths = isolines.collections[i].get_paths()
            for j in range(len(paths)):
                isoline_c1 = paths[j].vertices[:, 0] # shape: (*,)
                isoline_c2 = paths[j].vertices[:, 1] # shape: (*,)
                r = slice.convert2DCoordTo3DCoord(isoline_c1, isoline_c2) # shape: (*,3)
                x = r[:,0]
                y = r[:,1]
                z = r[:,2]
                level = levels[i]
                trace = go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(width=2, cmin=levels[0], cmax=levels[-1], color=np.ones_like(x)*level), showlegend=False)
                isolineDatList.append(trace)

                # annotation
                if cutIsolineNote(level):
                    continue
                c1i = isoline_c1[0]
                c2i = isoline_c2[0]
                c1f = isoline_c1[-1]
                c2f = isoline_c2[-1]
                length = np.sum(np.sqrt(np.diff(isoline_c1)**2 + np.diff(isoline_c2)**2))
                isLooped = (abs(c1i-c1f)<1e-10*sliceRepLength1) and (abs(c2i-c2f)<1e-10*sliceRepLength2)
                annotationList.append(dict(
                    text=level, x=x[0], y=y[0], z=z[0], font=dict(color='black'),
                    showarrow=bool((length<thresholdNoteArrow) and isLooped),
                    bgcolor='white', opacity=0.8
                ))

        return isolineDatList, annotationList

"""
上のクラスを統合して、cubeファイルのコンター図をプロットするための関数
"""
def visualizeCubeSlice(cube=None, cubeFile=None, outFile=None, slicePos=None, sliceNormal=None, pcaAuto=None, cmin=None, cmax=None, numIsoline=None, stepIsoline=None, cutIsolineNote=None, thresholdNoteArrow=None, cameraZoom=100, cameraRotate=0, printParam=False):
    import plotly.graph_objects as go

    # load cube data
    if cube is None:
        cube = Cube(filePath=cubeFile)
    unit = 'Bohr'

    node = cube.giveNodeCoord(unit=unit)
    xmin = node[:,:,:,0].min()
    xmax = node[:,:,:,0].max()
    ymin = node[:,:,:,1].min()
    ymax = node[:,:,:,1].max()
    zmin = node[:,:,:,2].min()
    zmax = node[:,:,:,2].max()

    cvis = CubeVisualizer()

    datList = []

    # set molecule plot
    if cube.giveMoleculeObj() is not None:
        molplotDatList = cvis.giveMoleculeSurfacePlot(cube, scale=1.5, unit=unit)
        datList.extend(molplotDatList)

    # set slice plot
    slice = Slice(cube, pos=slicePos, normal=sliceNormal, pcaAuto=pcaAuto, unit=unit)
    sliceDat = cvis.giveSlicePlot(slice, unit=unit)
    isolineDatList, annotationList = cvis.giveIsolinesPlot(slice, numIsoline=numIsoline, stepIsoline=stepIsoline, cutIsolineNote=cutIsolineNote, thresholdNoteArrow=thresholdNoteArrow, unit=unit)

    sliceDat['cmin'] = cmin
    sliceDat['cmax'] = cmax
    sliceDat['colorscale'] = 'rainbow'
    sliceDat['lighting'] = {'ambient':1.0}
    for d in isolineDatList:
        d['line']['colorscale'] = [[0,'rgb(0,0,0)'],[1,'rgb(0,0,0)']]

    datList.append(sliceDat)
    datList.extend(isolineDatList)

    # set camera parameter
    # convert to rad from deg
    cameraRotate = cameraRotate / 180 * np.pi
    sliceCenter = slice.giveSliceCenter(unit=unit)
    sliceNormal = slice.giveSliceNormalVector()
    sliceTangent = slice.giveSliceTangentVector()
    cameraPos = sliceNormal * (16.8 / np.abs(xmax-xmin) * 100 / cameraZoom)
    # rotation around normal axis
    cameraUp = sliceTangent * np.cos(cameraRotate) + sliceNormal * (sliceNormal@sliceTangent) * (1-np.cos(cameraRotate)) + np.cross(sliceNormal,sliceTangent)*np.sin(cameraRotate)

    fig = go.Figure(data=datList)
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=(xmin, xmax), showgrid=False, showticklabels=False, showbackground=False, title=''),
            yaxis=dict(range=(ymin, ymax), showgrid=False, showticklabels=False, showbackground=False, title=''),
            zaxis=dict(range=(zmin, zmax), showgrid=False, showticklabels=False, showbackground=False, title=''),
            aspectmode='manual',
            aspectratio=dict(x=1, y=(ymax-ymin)/(xmax-xmin), z=(zmax-zmin)/(xmax-xmin)),
            camera=dict(
                eye=dict(x=cameraPos[0], y=cameraPos[1], z=cameraPos[2]),
                up=dict(x=cameraUp[0], y=cameraUp[1], z=cameraUp[2])
            ),
            annotations=annotationList
        ),
        margin=dict(
            l=20, r=20, t=10, b=10
        )
    )

    if type(outFile) is str:
        fig.write_image(outFile)
    else:
        fig.show()

    if printParam:
        v = slice.giveSliceValue()

        print('cubeFile: {}'.format(cubeFile))
        print('minValue: {}'.format(np.nanmin(v)))
        print('maxValue: {}'.format(np.nanmax(v)))
        print('cmin: {}'.format(cmin))
        print('cmax: {}'.format(cmax))
        print('sliceCenter: {} {}'.format(sliceCenter,unit))
        print('sliceNormal: {}'.format(sliceNormal))
        print('sliceTangent: {}'.format(sliceTangent))
        print('cameraPos: {}'.format(cameraPos))
        print('cameraUp: {}'.format(cameraUp))

