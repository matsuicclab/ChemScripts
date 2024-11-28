import re

from rdkit import Chem

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
        symblist, _ = self.giveMolecularGeometry(orientation='Input')
        numAtom = len(symblist)
        return numAtom


    def giveMolecularGeometry(self, step=-1, unit='Angstrom', orientation=None):
        """
        step: step of geometry optimization to obtain (1,2,3,4,...,-1)
        unit: Angstrom or a.u.

        return: (symblist, xyzlist)
        """
        if unit == 'a.u.':
            factor = 1.889726126
        elif unit == 'Angstrom':
            factor = 1
        else:
            raise ValueError('Unsupported unit: {}'.format(unit))

        if step >= 1:
            step -= 1

        if orientation not in ['Input', 'Standard']:
            raise ValueError('Invalid orientation setting: {}'.format(orientation))

        idxlist_orient = [i for i,s in enumerate(self.__logdata) if '{} orientation:'.format(orientation) in s]
        if len(idxlist_orient) == 0:
            return None, None
        idx_orient = idxlist_orient[step]

        pedtab = Chem.GetPeriodicTable()
        symblist = []
        xyzlist = []
        countHyphenLine = 0
        for line in self.__logdata[idx_orient:]:
            if re.fullmatch('^-+$', line):
                countHyphenLine += 1
                continue
            if countHyphenLine == 3:
                break
            if countHyphenLine == 2:
                # center number, atomic number, atomic type, coord x, coord y, coord z (angstrom)
                _, n, _, x, y, z = line.split()

                # convert atomic number to element symbol
                s = pedtab.GetElementSymbol(n)
                # convert unit
                x *= factor
                y *= factor
                z *= factor

                # append
                symblist.append(s)
                xyzlist.append([float(x),float(y),float(z)])

        return symblist, xyzlist


    def giveXYZBlock(self, step=-1, unit='Angstrom', orientation=None):
        symblist, xyzlist = self.giveMolecularGeometry(step=step, unit=unit, orientation=orientation)
        xyzblock = '\n'.join(['{} {} {} {}'.format(s,x,y,z) for s, (x,y,z) in zip(symblist,xyzlist)])
        return xyzblock


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




