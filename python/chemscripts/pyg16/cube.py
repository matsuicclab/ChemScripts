import re
import itertools
import copy

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull

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


    def __init__fromFile(self, filePath=None, charge=None, comment=None):
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
        # (cubegenでの作成時にBohrとangstromを間接的に指定するがそれはcubegenの入力パラメータの単位を指している)
        # ここではunit='Bohr'として扱えばよく、一部を単位変換することはしない。
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

        self.__valueDim = valueDim

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
        self.__molecule = Molecule(atomicnumList=atomicNumData, xyzList=atomXYZData, unit='Bohr', charge=charge)

        # cubeデータ取り出し
        # 平坦化し、3次元の行列に変換
        self.__cubeData = np.array(
                            list(itertools.chain.from_iterable(data[6+numAtom:]))
                        ).reshape(*numGridPoint,valueDim)

        self.__sourceFilePath = filePath

        if comment is None:
            comment = 'load from {}'.format(filePath)
        self.setComment(comment)


    def __init__fromCubeData(self, cubeGrid=None, cubeData=None, moleculeObj=None, comment=None):
        """
        格子データが既に存在する場合に利用する

        cubeDataはnp.ndarrayであり、
        shapeは1次元(nx*ny*nz*valueDim,)または3次元(nx,ny,nz)(valueDim==1)、
        複数の格子データが存在する場合は4次元(nx,ny,nz,valueDim)とする
        """
        # 値チェック
        if cubeGrid is None or type(cubeGrid) is not CubeGrid:
            raise ValueError('type of cubeGrid must be chemscripts.pyg16.cube.CubeGrid')

        if cubeData is None:
            raise ValueError('cubeData is None')

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
        # shapeを(*numGridPoint, valueDim)に変形
        if len(cubeData.shape) == 1 and float(cubeData.shape[0]/np.prod(numGridPoint)).is_integer():
            cubeData = cubeData.reshape(*numGridPoint, -1)
        elif len(cubeData.shape) == 2 and float(cubeData.shape[0]/np.prod(numGridPoint)).is_integer():
            cubeData = cubeData.reshape(*numGridPoint, cubeData.shape[1])
        elif len(cubeData.shape) == 3 and cubeData.shape == tuple(numGridPoint):
            cubeData = cubeData.reshape(*numGridPoint, 1)
        elif len(cubeData.shape) == 4 and cubeData.shape[:-1] == tuple(numGridPoint):
            pass # reshapeの必要なし
        else:
            # いずれでもなかった場合
            raise ValueError('shape of cubeData must be (nx*ny*nz*valueDim,) or (nx,ny,nz)(valueDim==1), or (nx,ny,nz,valueDim)')

        # valueDim決定
        valueDim = int(np.prod(cubeData.shape) / np.prod(numGridPoint))

        self.__cubeGrid = cubeGrid
        self.__cubeData = cubeData
        self.__valueDim = valueDim

        # moleculeObj設定
        self.setMoleculeObj(moleculeObj)

        self.__sourceFilePath = None

        if comment is None:
            comment = 'generate from cubeData'
        self.setComment(comment)

    def __str__(self):
        return '\n'.join([
                    repr(self),
                    'comment: {}'.format(self.__comment),
                    'cubeGrid:', *[re.sub('^', '\t', s) for s in format(self.__cubeGrid).split('\n') if s != ''],
                    'moleculeObj:', *[re.sub('^', '\t', s) for s in format(self.__molecule).split('\n') if s != ''],
                    ''
                ])

    def __pos__(self):
        """
        + self
        """
        return Cube(cubeGrid=self.__cubeGrid, cubeData=self.__cubeData, moleculeObj=self.__molecule, comment=self.__comment)

    def __neg__(self):
        """
        - self
        """
        return Cube(cubeGrid=self.__cubeGrid, cubeData=-self.__cubeData, moleculeObj=self.__molecule, comment='Negation: (Cube: {})'.format(self.__comment))

    def __add_and_sub(self, other, sign1, sign2):
        """
        sign1 * self + sign2 * other
        """
        operation = 'addition' if sign1 * sign2 == 1 else 'subtraction'

        if type(other) is Cube:
            sign1str = ''
            sign2str = '+' if sign2 == 1 else '-'
            # どっちもCubeの場合は self + other かself - otherしかないためselfの符号は+
            newComment = '{}: {}(Cube: {}) {} (Cube: {})'.format(operation, sign1str, self.__comment, sign2str, other.__comment)

            if self.__valueDim != other.__valueDim:
                raise ValueError('Number of dimensions in cubeData does not match: {} & {}'.format(self.__valueDim, other.__valueDim))
            elif self.__cubeGrid == other.__cubeGrid:
                # 同じグリッド上でcubeデータを保持している場合
                return Cube(cubeGrid=self.__cubeGrid, cubeData=sign1 * self.__cubeData + sign2 * other.__cubeData, moleculeObj=self.__molecule, comment=newComment)
            else:
                # グリッドが異なる場合
                # selfのグリッド上でotherの値を補間して和を計算する
                interpolatedCubeData = other.interpolate(self.giveNodeCoord(unit='Bohr').reshape(-1,3), unit='Bohr')
                interpolatedCubeData = interpolatedCubeData.reshape(*self.giveNumGridPoint(), self.__valueDim)
                newCubeData = sign1 * self.__cubeData + sign2 * interpolatedCubeData
                return Cube(cubeGrid=self.__cubeGrid, cubeData=newCubeData, moleculeObj=self.__molecule, comment=newComment)

        elif type(other) in [float, int, np.float16, np.float32, np.float64, np.int16, np.int32, np.int64]:
            sign1str = '' if sign1 == 1 else '-'
            sign2str = '+' if sign2 == 1 else '-'
            newComment = '{}: {}(Cube: {}) {} ({})'.format(operation, sign1str, self.__comment, sign2str, other)

            return Cube(cubeGrid=self.__cubeGrid, cubeData=sign1 * self.__cubeData + sign2 * other, moleculeObj=self.__molecule, comment=newComment)
        else:
            raise TypeError('invalid type: {}'.format(type(other)))

    def __add__(self, other):
        """
        self + other
        """
        return self.__add_and_sub(other, 1, 1)

    def __radd__(self, other):
        """
        other + self
        """
        return self.__add_and_sub(other, 1, 1)

    def __sub__(self, other):
        """
        self - other
        """
        return self.__add_and_sub(other, 1, -1)

    def __rsub__(self, other):
        """
        other - self
        """
        return self.__add_and_sub(other, -1, 1)

    def __mul_and_div(self, other, sign1, sign2):
        """
        self ^ sign1 * other ^ sign2
        """
        operation = 'multiply' if sign1 * sign2 == 1 else 'division'

        if type(other) is Cube:
            sign1str = ''
            sign2str = '*' if sign2 == 1 else '/'
            # どっちもCubeの場合は self * other かself / otherしかないためselfの符号は''
            newComment = '{}: {}(Cube: {}) {} (Cube: {})'.format(operation, sign1str, self.__comment, sign2str, other.__comment)

            if self.__valueDim != other.__valueDim:
                raise ValueError('Number of dimensions in cubeData does not match: {} & {}'.format(self.__valueDim, other.__valueDim))
            elif self.__cubeGrid == other.__cubeGrid:
                # 同じグリッド上でcubeデータを保持している場合
                return Cube(cubeGrid=self.__cubeGrid, cubeData=self.__cubeData**sign1 * other.__cubeData**sign2, moleculeObj=self.__molecule, comment=newComment)
            else:
                # グリッドが異なる場合
                # selfのグリッド上でotherの値を補間して和を計算する
                interpolatedCubeData = other.interpolate(self.giveNodeCoord(unit='Bohr').reshape(-1,3), unit='Bohr')
                interpolatedCubeData = interpolatedCubeData.reshape(*self.giveNumGridPoint(), self.__valueDim)
                newCubeData = self.__cubeData**sign1 * interpolatedCubeData**sign2
                return Cube(cubeGrid=self.__cubeGrid, cubeData=newCubeData, moleculeObj=self.__molecule, comment=newComment)

        elif type(other) in [float, int, np.float16, np.float32, np.float64, np.int16, np.int32, np.int64]:
            sign1str = '' if sign1 == 1 else '1/'
            sign2str = '*' if sign2 == 1 else '/'
            newComment = '{}: {}(Cube: {}) {} ({})'.format(operation, sign1str, self.__comment, sign2str, other)

            return Cube(cubeGrid=self.__cubeGrid, cubeData=self.__cubeData**sign1 * other**sign2, moleculeObj=self.__molecule, comment=newComment)
        else:
            raise TypeError('invalid type: {}'.format(type(other)))

    def __mul__(self, other):
        """
        self * other
        """
        return self.__mul_and_div(other, 1, 1)

    def __rmul__(self, other):
        """
        other * self
        """
        return self.__mul_and_div(other, 1, 1)

    def __truediv__(self, other):
        """
        self / other
        """
        return self.__mul_and_div(other, 1, -1)

    def __rtruediv__(self, other):
        """
        other / self
        """
        return self.__mul_and_div(other, -1, 1)

    def __pow__(self, other):
        """
        self ** other
        """
        if type(other) not in [int, float]:
            raise ValueError()
        newCubeData = self.__cubeData ** other
        newComment = '({})**({})'.format(self.__comment, other)
        return Cube(cubeGrid=self.__cubeGrid, cubeData=newCubeData, moleculeObj=self.__molecule, comment=newComment)

    def setComment(self, comment):
        if comment is None:
            comment = '' # 空文字にする
        elif type(comment) is not str:
            raise TypeError('type of comment must be str: {}'.format(type(comment)))
        elif '\n' in comment:
            raise ValueError('Multi-line string is invalid as comment')
        self.__comment = comment

    def setMoleculeObj(self, moleculeObj):
        if moleculeObj is not None and type(moleculeObj) is not Molecule:
            # Noneは許容するが型チェックだけしとく
            raise TypeError('type of moleculeObj must be chemscripts.molecule.Molecule')

        self.__molecule = moleculeObj

    def giveComment(self):
        return self.__comment

    def giveSourceFilePath(self):
        return self.__sourceFilePath

    def giveCubeData(self):
        """
        cubeデータのコピーを返す
        """
        return copy.deepcopy(self.__cubeData)

    def giveCubeGrid(self):
        return self.__cubeGrid

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
            numAtom = 0 if self.__molecule is None else self.__molecule.giveNumAtom()

            header1 = 'Cube Data generated by chemscripts.pyg16.cube.Cube.write()\n'
            header2 = '{}\n'.format(self.__comment)
            header3 = '{} {} {} {} {}\n'.format(numAtom, *self.giveStartingPoint(unit='Bohr'), self.__valueDim)
            header456 = ''.join(['{} {} {} {}\n'.format(n,*v) for n,v in zip(self.giveNumGridPoint(),self.giveStepVector(unit='Bohr'))])
            if self.__molecule is not None:
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
        elif 'coords' in args.keys():
            self.__init__fromCubicCoords(**args)
        elif 'startingPoint' in args.keys():
            self.__init__fromParam(**args)
        elif 'cubeGrid' in args.keys():
            self.__init__fromCubeGrid(**args)
        else:
            raise ValueError('args must contain moleculeObj or startingPoint or cubeGrid')

    def __init__fromMoleculeObj(self, moleculeObj=None, axesMethod=None, step=0.5, padding=3.0, unit='Angstrom'):
        """
        unit: unit of step and padding
        axesMethod: 'Direct' or 'PCA' or 'PCA-ignoreHs'
        """

        # molecule
        if moleculeObj is None:
            raise ValueError()
        if type(moleculeObj) is not Molecule:
            raise TypeError()
        # step
        if type(step) not in [int, float, np.float64, np.float32]:
            raise TypeError('step must be a number: type is {}'.format(type(step)))
        # padding
        if type(padding) not in [int, float, np.float64, np.float32]:
            raise TypeError('padding must be a number: type is {}'.format(type(padding)))
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
        # それぞれのベクトルの長さをstepに調整
        v1 = v1 / np.linalg.norm(v1) * step
        v2 = v2 / np.linalg.norm(v3) * step
        v3 = v3 / np.linalg.norm(v3) * step
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
        startingPoint = (V @ minijkAtom.T).T - padding * (v1 + v2 + v3) / step
        endingPoint = (V @ maxijkAtom.T).T + padding * (v1 + v2 + v3) / step

        self.__init__fromParam(startingPoint=startingPoint, stepVector=stepVector, endingPoint=endingPoint, unit=unit)


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

    def __init__fromCubicCoords(self, coords=None, unit=None):
        # 最初に凸包を考えて全体の平行六面体の頂点座標を割り出す
        hull = ConvexHull(coords)
        vertices = coords[hull.vertices] # shape: (6,3)
    
        # 全体平行六面体の最も離れた頂点同士をstartingPointとendingPointにする
        # 幾つか候補がある場合にはendingPoint-startingPointがベクトル[1,1,1]に近い向きになるペアを使用
        dist2 = cdist(vertices, vertices)
        distmax_pairidxs = np.where(dist2 == dist2.reshape(-1)[np.argmax(dist2)]) # shape: (numCandidates, 2)
        vs = np.array([vertices[i]-vertices[j] for i,j in zip(*distmax_pairidxs)]) # shape: (numCandidates, 3)
        # 候補から一つだけ選び出す
        overlaps = vs / np.linalg.norm(vs, axis=1)[:,np.newaxis] @ np.array([1,1,1]) # shape: (numCandidates,)
        endingVertexIdx   = distmax_pairidxs[0][np.argmax(overlaps)]
        startingVertexIdx = distmax_pairidxs[1][np.argmax(overlaps)]
        endingPoint = vertices[endingVertexIdx]
        startingPoint = vertices[startingVertexIdx]
        v = endingPoint - startingPoint
    
        # 全体平行六面体の各頂点とstartingPoint間の距離を計算
        dist = np.linalg.norm(vertices - startingPoint, axis=1)
        # 1番目と2番目に近い距離にある頂点をstartingPointと隣り合う頂点として選定
        # (三角不等式により2番目までは隣り合う頂点の選定をすることができる)
        v1 = vertices[np.argsort(dist)[1]] - startingPoint
        v2 = vertices[np.argsort(dist)[2]] - startingPoint
        v3 = v - v1 - v2
    
        # v1,v2,v3の順番を調整する
        #  - v1ができるだけx軸方向になるようにする
        #  - [v1,v2,v3]の行列式が正になるようにする (右手系)
        V = np.array([v1,v2,v3])
        V = V[np.argsort(-1 * V / np.linalg.norm(V, axis=1)[:,np.newaxis] @ np.array([1,0,0]))]
        if np.linalg.det(V) < 0:
            V = V[[0,2,1]]
        v1, v2, v3 = V
    
        # 各軸方向のグリッド数を推定する
        # 各グリッドの座標を[v1,v2,v3]で展開すると座標成分は[i/(N1-1), j/(N2-1), k/(N3-1)] (i,j,k in [0,N-1])
        # この2乗和からグリッド数を推定することができる。
        c = np.linalg.solve(V.T, (coords-startingPoint).T)
        N = coords.shape[0]
        numGridPoint = 1 + 1 / (6 / N * np.sum(c**2, axis=1) - 2)
        numGridPoint = np.round(numGridPoint).astype(np.int32)
        if np.prod(numGridPoint) != N:
            raise ValueError('Invalid numGridPoint: {} != {} * {} * {}'.format(N, *numGridPoint))
        
        # メンバ変数に設定
        self.__init__fromParam(startingPoint=startingPoint, numGridPoint=numGridPoint, endingPoint=endingPoint, unit=unit)
        

    def __init__fromCubeGrid(self, cubeGrid=None, stepDetailRatio=1, numMarginGrid=0, shiftGrid=0):
        """
        stepDetailRatio: 非ゼロ整数
                    newStepVector = oldStepVector / (abs(ratio) ** sign(ratio))
        numMarginGrid: 整数
                    与えられたcubeGridの外側にどれだけ広げるか(newStepVector単位)
        shiftGrid: [0,1)の実数
                    startingPointをどれだけズラすか(newStepVector単位)
        """
        if type(cubeGrid) is not CubeGrid:
            raise TypeError('type of cubeGrid must be CubeGrid')
        if type(stepDetailRatio) is not int:
            raise TypeError('type of stepDetailRatio must be int')
        if type(numMarginGrid) is not int:
            raise TypeError('type of numMarginGrid must be int')
        if type(shiftGrid) not in [int, float]:
            raise TypeError('shiftGrid must be a number')

        if stepDetailRatio == 0:
            raise ValueError('stepDetailRatio must not be zero')
        if shiftGrid < 0 or shiftGrid >= 1:
            raise ValueError('shiftGrid must be a real number greater than or equal to 0 and less than 1')

        # 単位決定
        newUnit = cubeGrid.__unit

        # stepVector決定
        # numGridPoint決定
        if stepDetailRatio > 0:
            newStepVector = cubeGrid.__stepVector / stepDetailRatio
            newNumGridPoint = (stepDetailRatio * (cubeGrid.__numGridPoint - 1) + 1) + numMarginGrid * 2
        else:
            newStepVector = cubeGrid.__stepVector * np.abs(stepDetailRatio)
            newNumGridPoint = ((cubeGrid.__numGridPoint - 1) // np.abs(stepDetailRatio) + 1) + numMarginGrid * 2
        v1, v2, v3 = newStepVector

        # startingPoint決定
        newStartingPoint = cubeGrid.__startingPoint + (v1+v2+v3) * (- numMarginGrid + shiftGrid)

        # 構築
        self.__init__fromParam(startingPoint=newStartingPoint, stepVector=newStepVector, numGridPoint=newNumGridPoint, unit=newUnit)


    def __str__(self):
        return '\n'.join([
                    repr(self),
                    'startingPoint: {}'.format(self.__startingPoint),
                    'endingPoint: {}'.format(self.__endingPoint),
                    'stepVector: {}'.format(self.__stepVector),
                    'numGridPoint: {}'.format(self.__numGridPoint),
                    'unit: {}'.format(self.__unit),
                    ''
                ])

    def __eq__(self, other):
        """
        self == other
        """
        if type(other) is not CubeGrid:
            # 型判定
            return False
        elif self is other:
            # 同値判定
            return True
        elif (self.__startingPoint == other.__startingPoint).all() and \
            (self.__endingPoint == other.__endingPoint).all() and \
            (self.__stepVector == other.__stepVector).all() and \
            (self.__numGridPoint == other.__numGridPoint).all() and \
            self.__unit == other.__unit:
            # 等価判定
            return True
        else:
            return False

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
        unit: unit of r
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
        unit: unit of r (return)
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

    def giveDeltaV(self, unit=None):
        """
        グリッド一間隔辺りの体積を計算
        """
        return np.abs(np.linalg.det(self.giveStepVector(unit=unit)))


"""
上のクラスを統合して、cubeファイルのコンター図をプロットするための関数

この関数は以前はPlotlyを用いて実装されていたが、
現在はPyVistaを用いた実装に移されている。
Plotlyでスライスや等値線を実装することはできていたが、
コードの量が非常に多く保守しにくいため、PyVistaに移行している。
移行後の関数はchemscripts.pyg16.plotcube.visualizeCubeSlice()である。
"""
def visualizeCubeSlice(cube=None, cubeFile=None, outFile=None, slicePos=None, sliceNormal=None, pcaAuto=None, cmin=None, cmax=None, numIsoline=None, stepIsoline=None, cutIsolineNote=None, noteFormat=None, showNote=True, thresholdNoteArrow=None, cameraZoom=100, cameraRotate=0, printParam=False):
    raise NotImplementedError()


