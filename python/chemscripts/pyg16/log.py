import re

from rdkit import Chem

from chemscripts.molecule import Molecule
from chemscripts.unit import getUnitConversionFactor

class Log:
    def __init__(self, **args):
        if 'filePath' in args.keys():
            self.__init__fromFile(**args)
        elif 'logData' in args.keys():
            self.__init__fromLogData(**args)
        else:
            raise ValueError('args must contain filePath or cubeData')

    def __init__fromFile(self, filePath=None):
        """
        load log file
        """
        # ファイル読み込み
        with open(filePath, mode='r') as f:
            self.__logdata = [s.strip('\n') for s in f.readlines()]

    def __init__fromLogData(self, logData=None):
        self.__logdata = logData


    def giveNumAtom(self):
        """
        return: the number of atoms (int)
        """
        idx = [i for i,s in enumerate(self.__logdata) if 'NAtoms= ' in s][0]
        numAtom = re.sub(' .+$', '', re.sub('^[^=]+= +', '', self.__logdata[idx]))
        return int(numAtom)


    def giveMoleculeObj(self, step=-1, orientation=None):
        """
        step: step of geometry optimization to obtain (0,1,2,3,...,-1)

        return: chemscripts.molecule.Molecule
        """
        xyzBlock = self.giveXYZBlock(step=step, unit='Angstrom', orientation=orientation)
        if xyzBlock is None:
            return None
        molecule = Molecule(xyzBlock=xyzBlock, unit='Angstrom')

        return molecule


    def giveXYZBlock(self, step=-1, unit='Angstrom', orientation=None):
        """
        step: step of geometry optimization to obtain (1,2,3,4,...,-1)
        unit: Angstrom or a.u.

        return: string of xyzBlock
        """
        factor = getUnitConversionFactor(oldunit='Angstrom', newunit=unit)

        if type(step) is not int:
            raise TypeError('type of step must be int')

        if orientation not in ['Input', 'Standard']:
            raise ValueError('Invalid orientation setting: {}'.format(orientation))

        idxlist_orient = [i for i,s in enumerate(self.__logdata) if '{} orientation:'.format(orientation) in s]
        if len(idxlist_orient) == 0:
            return None
        idx_orient = idxlist_orient[step]

        table = Chem.GetPeriodicTable()
        xyzBlock = []
        countHyphenLine = 0
        for line in self.__logdata[idx_orient:]:
            if re.fullmatch('^ -+$', line):
                countHyphenLine += 1
                continue
            if countHyphenLine == 3:
                break
            if countHyphenLine == 2:
                # center number, atomic number, atomic type, coord x, coord y, coord z (angstrom)
                _, n, _, x, y, z = line.split()

                # convert atomic number to element symbol
                s = table.GetElementSymbol(int(n))
                # convert unit
                x = float(x) * factor
                y = float(y) * factor
                z = float(z) * factor

                # append
                xyzBlock.append('{} {} {} {}'.format(s,x,y,z))

        xyzBlock = '\n'.join([str(len(xyzBlock)), 'comment'] + xyzBlock)
        return xyzBlock


    def givePCMModel(self):
        idxlist = [i for i,s in enumerate(self.__logdata) if 'Polarizable Continuum Model (PCM)' in s]
        if len(idxlist) == 0:
            return None
        return re.sub('\.$', '', re.sub('^.+: ','', self.__logdata[idxlist[-1]+2]))

    def givePCMAtomicRadii(self):
        idxlist = [i for i,s in enumerate(self.__logdata) if 'Polarizable Continuum Model (PCM)' in s]
        if len(idxlist) == 0:
            return None
        return re.sub('\.$', '', re.sub('^.+: ','', self.__logdata[idxlist[-1]+3]))

    def givePCMCavityType(self):
        idxlist = [i for i,s in enumerate(self.__logdata) if 'Polarizable Continuum Model (PCM)' in s]
        if len(idxlist) == 0:
            return None
        return re.sub(' +(.+$', '', re.sub('^.+: ','', self.__logdata[idxlist[-1]+7]))

    def givePCMCavityScalingFactor(self):
        idxlist = [i for i,s in enumerate(self.__logdata) if 'Polarizable Continuum Model (PCM)' in s]
        if len(idxlist) == 0:
            return None
        return re.sub(')\.$', '', re.sub('^.+([^=]+=','', self.__logdata[idxlist[-1]+7]))

    def givePCMEpsilon(self):
        idxlist = [i for i,s in enumerate(self.__logdata) if 'Polarizable Continuum Model (PCM)' in s]
        if len(idxlist) == 0:
            return None
        return re.sub(' .+$', '', re.sub('^.+Eps= *','', self.__logdata[idxlist[-1]+18]))




