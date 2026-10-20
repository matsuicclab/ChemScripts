import sys
import re

import numpy as np
import pyvista as pv

from chemscripts.pyg16.cube import Cube
from chemscripts.unit import checkInvalidUnit, getUnitConversionFactor



def generate_cubemesh(cube, unit=None):
    """
    cubeオブジェクトからPyVistaのmeshオブジェクトを構築
    cubeDataはpoint_data['cubedata']に格納
    
    args:
        cube:
            chemscripts.pyg16.cube.Cube
            cubeオブジェクト
        
    return:
        cubemesh:
            pyvista.core.pointset.StructuredGrid
            cubeオブジェクトを変換して得られるmeshオブジェクト
    """
    # cubeデータをmeshに変換
    startingPoint = cube.giveStartingPoint(unit=unit)
    stepVector = cube.giveStepVector(unit=unit)
    numGridPoint = cube.giveNumGridPoint()

    v1, v2, v3 = stepVector
    n1, n2, n3 = numGridPoint
    I, J, K = np.meshgrid(np.arange(n1), np.arange(n2), np.arange(n3), indexing='ij')
    R = startingPoint + I[:,:,:,np.newaxis] * v1 + J[:,:,:,np.newaxis] * v2 + K[:,:,:,np.newaxis] * v3
    X, Y, Z = R[:,:,:,0], R[:,:,:,1], R[:,:,:,2]
    cubemesh = pv.StructuredGrid(X, Y, Z)
    cubemesh.point_data['cubedata'] = cube.giveCubeData().reshape(cubemesh.n_points, -1)

    return cubemesh

def generate_slicemesh(cube, pos=None, normal=None, pcaAuto=None, unit=None):
    """
    cubeオブジェクトから断面(スライス)をPyVistaのmeshオブジェクトとして構築
    
    args:
        cube:
            chemscripts.pyg16.cube.Cube
            cubeオブジェクト
        pos:
            np.ndarray (shape: (3,))
            スライスの位置ベクトル
        normal:
            np.ndarray (shape: (3,))
            スライスの法線ベクトル
        pcaAuto:
            bool
            分子の原子座標についてPCAを実行し、その第1軸・第2軸でスライス面を定義する場合はTrue
        unit:
            str
            posの単位系 ('Angstrom' / 'Bohr')
    
    return:
        slicemesh:
            pyvista.core.pointset.PolyData
            スライスを表すPyVistaのmeshオブジェクト
        pos:
            np.ndarray (shape: (3,))
            スライスの位置ベクトル
        normal:
            np.ndarray (shape: (3,))
            スライスの法線ベクトル
    """
    # cubicなmeshを生成
    cubemesh = generate_cubemesh(cube, unit=unit)

    if pcaAuto:
        # 原子核の配置からスライス面を決定
        molecule = cube.giveMoleculeObj()
        # TODO: PCA-ignoreHsだと期待通りの面にならない場合もあるので要再考
        pos, _, _, normal = molecule.generateStandardizedCoordSystem(unit=unit,method='PCA')

    # 型チェック
    if pos is None or normal is None:
        raise ValueError()
    elif type(pos) is not np.ndarray or type(normal) is not np.ndarray:
        raise TypeError()
    elif pos.shape != (3,) or normal.shape != (3,):
        raise ValueError() 

    # スライス面meshを生成
    slicemesh = cubemesh.slice(normal=normal, origin=pos)
    return slicemesh, pos, normal

def generate_isolines(slicemesh, numIsoline=None, stepIsoline=1, labelFormat=None, cutIsolineLabel=None, epsz=1e-2):
    """
    スライス上の等値線をPyVistaのmeshオブジェクトとして生成
    
    args:
        numIsoline:
            int
            等値線の数を設定
        stepIsoline:
            float
            等値線の間隔を設定
        labelFormat:
            str
            labelの数値のformatを設定
        cutIsolineLabel:
            function
            labelを非表示にしたいlevelを関数で指定
        epsz:
            float
            等値線とスライスの間の距離
            epsz=0でスライス上に等値線を配置することになるが、重なってしまうために正常に描画されないため、基本epsz!=0とする
        
    return:
        frontIsolineList:
            list of pyvista.core.pointset.PolyData 
            スライスの表面に配置する等値線のリスト
        backIsolineList:
            list of pyvista.core.pointset.PolyData
            スライスの裏面に配置する等値線のリスト
        labelDictList:
            list of dict
            等値線の値をラベルとして表示する際に使う設定のリスト
            pyvista.Plotter.add_point_labels()に入力する
    """
    if cutIsolineLabel is None:
        cutIsolineLabel = lambda level: False
    if labelFormat is None:
        labelFormat = ''
    labelFormat = '{' + labelFormat + '}'

    # isolineのレベルを決定
    level_min = np.nanmin(slicemesh.point_data['cubedata'])
    level_max = np.nanmax(slicemesh.point_data['cubedata'])
    if numIsoline is not None:
        # numIsolineが指定されている場合はレベルの数からレベルを決定する
        levels = np.linspace(level_min, level_max, numIsoline)
    else:
        # stepIsolineが指定されている場合はレベルの間隔を指定する
        # levelがstepIsolineの整数倍になるように調整
        maxNumIsoline = 100
        if (level_max-level_min)/stepIsoline > maxNumIsoline:
            if np.abs(level_max) > np.abs(level_min):
                adjusted_level_max = level_min + maxNumIsoline * stepIsoline
                adjusted_level_min = level_min
            else:
                adjusted_level_max = level_max
                adjusted_level_min = level_max - maxNumIsoline * stepIsoline
        else:
            adjusted_level_max = level_max
            adjusted_level_min = level_min

        levels = np.arange(
                    np.ceil(adjusted_level_min/stepIsoline)*stepIsoline,
                    (np.floor(adjusted_level_max/stepIsoline)+1)*stepIsoline,
                    stepIsoline
                )

    # sliceの面の法線を計算しておく -> 'Normals'
    # これを使ってisolineをsliceから少しずらした位置に配置するようにする
    slicemesh.compute_normals(cell_normals=False, point_normals=True, inplace=True)

    frontIsolineList = []
    backIsolineList = []
    labelDictList = []
    for level in levels:
        isoline = slicemesh.contour(isosurfaces=[level])
        if isoline.n_points == 0:
            continue

        # 等値線の色は黒で固定する
        isoline.point_data['rgb'] = np.zeros((isoline.n_points, 3))

        # isolineListに追加する前に
        # slice面の表方向と裏方向にすこしだけズラす
        # これをしないとsliceとisolineが重なっているために綺麗に線が描画されない
        normal = isoline.point_data['Normals']
        isoline1 = isoline.copy()
        isoline2 = isoline.copy()
        isoline1.points += epsz * normal
        isoline2.points -= epsz * normal
        frontIsolineList.append(isoline1)
        backIsolineList.append(isoline2)

        if cutIsolineLabel(level):
            # ラベルを非表示にしたいlevelであればラベルの追加処理を行わない
            continue
        # 等値線のラベルを追加
        labelcoord = isoline.points[len(isoline.points)//2]
        labelDictList.append({
            "points": labelcoord[np.newaxis,:],
            "labels": [labelFormat.format(level)],
            "point_size": 0,
            "shadow": True,
            "shape": None,
            "shape_color": "white",
            "justification_horizontal": 'center',
            "justification_vertical": 'center',
            # 本当は"shape": "rect"で文字の後ろに白い四角を配置したかったが、
            # 文字の中心が指定した座標にくるように(justification: center)設定すると、
            # 四角がズレてしまう(PyVistaのバグ)ので、
            # 一旦文字に影をつけて対応する
        })

    return frontIsolineList, backIsolineList, labelDictList

def atomicnum_to_radius(n, scale=1.0, unit=None):
    """
    原子番号リストから原子半径リストを与える
    
    原子半径はGaussViewのデフォルトを参考に設定:
        周期      2r[A]
        1(1-2)       0.7878
        2(3-10)      1.0918
        3(11-18)     1.1801
        4(19-36)     1.2477
        5(37-54)     1.29
        6(55-86)     1.34
        7(87-118)    1.3903

    args:
        n:
            list of int
            原子番号リスト
        scale:
            float
            原子半径のスケール
            半径を大きくしたい場合は1.0より大きい数を指定する
        unit:
            str
            戻り値の原子半径の単位の指定

    return:
        radius:
            np.ndarray (shape: (numAtom,))
            各原子の原子半径リスト
    """
    def __func(n):
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

    return np.array([__func(ni) for ni in n]) * scale * factor

def atomicnum_to_rgb(n):
    """
    原子番号から原子の色を設定
    色はGaussViewのデフォルトを参考に設定
    
    args:
        n:
            list of int
            原子番号リスト
    return:
        rgb:
            np.ndarray (shape: (numAtom, 3))
            原子の色のリスト
            r,g,bの順であり、値は0 ~ 1に規格化されている
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
    rgb = [colorDict[ni] for ni in n] # shape: (numAtom, 3)

    return np.array(rgb) / 255

def generate_molmesh(molecule, addBonds=False, bondScale=1.0, unit=None):
    """
    moleculeオブジェクトからPyVistaのmeshオブジェクトを構築する
    
    args:
        molecule:
            chemscripts.molecule.Molecule
            変換するmoleculeオブジェクト
        addBonds:
            bool
            結合を表現するmeshオブジェクトを構築する(True)かしない(False)かを指定
        bondScale:
            float
            結合を表現する円柱の半径のscaleを指定
        unit:
            str
            meshオブジェクトのpoint座標の単位を指定
    
    return:
        atommeshList:
            list of pyvista.core.pointset.PolyData
            原子を表すPyVistaの球meshオブジェクトのリスト
        bondmeshList:
            list of pyvista.core.pointset.PolyData
            結合を表すPyVistaの円柱meshオブジェクトのリスト
    """
    # 分子構造の可視化
    atomicnumList = molecule.giveAtomicnumList()
    centerList = molecule.giveXYZArray(unit=unit)
    radiusList = atomicnum_to_radius(atomicnumList, unit=unit)
    rgbList = atomicnum_to_rgb(atomicnumList)

    # 原子核を追加
    atommeshList = []
    for center, radius, rgb in zip(centerList, radiusList, rgbList):
        sphere = pv.Sphere(radius=radius, center=center)
        sphere.point_data['rgb'] = rgb * np.ones((sphere.n_points, 3))
        atommeshList.append(sphere)

    # 結合を追加
    bondmeshList = []
    if addBonds:
        # 結合の棒の半径
        w = 0.05 * getUnitConversionFactor('Angstrom', unit) * bondScale
        # 結合の色
        rgb = np.array([150,150,150]) / 255

        # rdkitのMolオブジェクトを取得し、結合情報を復元
        mol = molecule.generateRDKitMolObj()
        bonds = mol.GetBonds()
        for b in bonds:
            atomIdx1 = b.GetBeginAtomIdx()
            atomIdx2 = b.GetEndAtomIdx()
            bondtype = str(b.GetBondType())

            r1 = centerList[atomIdx1]
            r2 = centerList[atomIdx2]

            bondCenter = (r1 + r2) / 2
            bondDirection = r2 - r1
            bondHeight = np.linalg.norm(bondDirection)
            cylinder = pv.Cylinder(center=bondCenter, direction=bondDirection, radius=w, height=bondHeight)
            cylinder.point_data['rgb'] = rgb * np.ones((cylinder.n_points, 3))
            bondmeshList.append(cylinder)

    return atommeshList, bondmeshList


def visualizeCubeSlice(
        cube=None, cubeFile=None, outFile=None, 
        slicePos=None, sliceNormal=None, pcaAuto=None, 
        cmin=None, cmax=None, 
        numIsoline=None, stepIsoline=None, cutIsolineLabel=None, 
        labelFormat=None, showLabels=True,
        printParam=False
    ):
    """
    cubeデータを入力し、
    スライス上のコンター(と等値線)と分子構造を可視化する
    
    args:
        cube:
            chemscripts.pyg16.cube.Cube
            可視化するcubeオブジェクト
        cubeFile:
            str
            可視化するcubeファイルへのパス
        outFile:
            str
            可視化結果の出力先パス
            Noneを指定した場合、GUIウインドウ上に描画を試みる
        slicePos:
            np.ndarray (shape: (3,))
            スライスの位置ベクトル
        sliceNormal:
            np.ndarray (shape: (3,))
            スライスの法線ベクトル
        pcaAuto:
            bool
            スライスの位置と法線を分子構造のPCA結果をもとに決定するかどうか
        cmin:
            float
            カラーバーの最小値
        cmax:
            float
            カラーバーの最大値
        numIsoline:
            int
            等値線の数
        stepIsoline:
            float
            等値線のレベル間隔
        cutIsolineLabel:
            function
            labelを非表示にしたいlevelを関数で指定
            lambda level: (level > 10) | (level < -10)
            など
        labelFormat:
            str
            labelの数値のformat指定: ':.2f'など
        showLabels:
            bool
            等値線のレベルのラベルを表示するかどうか
        printParam:
            bool
            デバッグ用
            
    """
    # load cube data
    if cube is None:
        cube = Cube(filePath=cubeFile)
    else:
        cubeFile = cube.giveSourceFilePath()
    
    unit = 'Bohr'

    # スライス
    slicemesh, slicePos, sliceNormal = generate_slicemesh(cube=cube, pos=slicePos, normal=sliceNormal, pcaAuto=pcaAuto, unit=unit)

    # 等値線
    frontIsolineList, backIsolineList, labelDictList = generate_isolines(slicemesh, numIsoline=numIsoline, stepIsoline=stepIsoline, labelFormat=labelFormat, cutIsolineLabel=cutIsolineLabel)

    # 分子構造
    molecule = cube.giveMoleculeObj()
    if molecule is not None:
        atommeshList, bondmeshList = generate_molmesh(molecule, unit=unit, addBonds=True)
    else:
        atommeshList, bondmeshList = [], []

    # Plotter用意
    window_size = (1440, 810)
    if outFile is None:
        # 出力先が指定されていない場合はGUIウインドウに描画を試みる
        pl = pv.Plotter(window_size=window_size)
    else:
        # 出力先が指定されている場合は画像出力を試みる
        pl = pv.Plotter(off_screen=True, window_size=window_size)
    pl.add_axes()

    # cubeのスライスと等値線を描画
    pl.add_mesh(
        slicemesh, 
        scalars='cubedata', cmap='turbo', clim=[cmin, cmax], 
        scalar_bar_args={'title': ''},
        lighting=False
    )
    for frontIsoline, backIsoline, labelDict in zip(frontIsolineList, backIsolineList, labelDictList):
        pl.add_mesh(frontIsoline, scalars='rgb', rgb=True, line_width=2)
        pl.add_mesh(backIsoline, scalars='rgb', rgb=True, line_width=2)
        if showLabels:
            pl.add_point_labels(**labelDict)
    # 原子の球と結合の円柱を描画
    # これらはローポリで綺麗に表示したいので、スムーズシェードを有効にする
    for atommesh in atommeshList:
        pl.add_mesh(atommesh, scalars='rgb', rgb=True, smooth_shading=True, ambient=0.3)
    for bondmesh in bondmeshList:
        pl.add_mesh(bondmesh, scalars='rgb', rgb=True, smooth_shading=True, ambient=0.3)

    # カメラの位置設定
    focal_point = slicePos
    # カメラとmeshの間の距離を決定: スライスの幅に応じて距離をとるように設定
    distance = 4.0 * np.max(np.linalg.norm(slicemesh.points - focal_point, axis=1))
    camera_position = slicePos + distance * sliceNormal
    viewup = np.array([0,0,1])
    viewup = viewup - (viewup @ sliceNormal) * sliceNormal

    viewup = viewup / np.linalg.norm(viewup)
    #TODO: 回転角での指定の方がよくない？
    cp = pv.CameraPosition(position=camera_position, focal_point=focal_point, viewup=viewup)
    pl.camera_position = cp
    pl.enable_parallel_projection() # 平行投影: この呼び出しはカメラの位置を指定した後でなければ正しく反映されない

    if outFile is None:
        # GUIウインドウへ
        pl.show()
    else:
        # 画像ファイルへ
        pl.screenshot(outFile)
        
    if printParam:
        print('cubeFile: {}'.format(cubeFile))
        level_min = np.nanmin(slicemesh.point_data['cubedata'])
        level_max = np.nanmax(slicemesh.point_data['cubedata'])
        print('minValue: {}'.format(level_min))
        print('maxValue: {}'.format(level_max))
        print('cmin: {}'.format(cmin))
        print('cmax: {}'.format(cmax))
        print('slice.position: {} {}'.format(slicePos,unit))
        print('slice.normal: {}'.format(sliceNormal))
        print('camera.position: {}'.format(camera_position))
        print('camera.viewup: {}'.format(viewup))
        

if __name__ == '__main__':
    cubeFile = sys.argv[1]
    outFile = re.sub('\.cube?$', '', cubeFile) + '.png'
    visualizeCubeSlice(
        cubeFile=cubeFile,
        outFile=outFile,
        cmin=None, cmax=None,
        stepIsoline=0.5,
        pcaAuto=True,
        labelFormat=':.2f'
    )
    
    
    