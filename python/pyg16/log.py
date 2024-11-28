import re

from rdkit import Chem

class Log:
    def __init__(self, filePath):
        """
        load log file
        """
        # ファイル読み込み
        with open(filePath, mode='r') as f:
            self.__logdata = [s.strip('\n') for s in f.readlines()]

    def giveNumAtom(self):
        symblist, _ = self.giveMolecularGeometry_in_InputOrientation()
        numAtom = len(symblist)
        return numAtom

    def giveMolecularGeometry_in_InputOrientation(self, step=-1, unit='Angstrom'):
        """
        step: step of geometry optimization to obtain (1,2,3,4,...,-1)
        unit: Angstrom or a.u.
        """
        if unit == 'a.u.':
            factor = 1.889726126
        elif unit == 'Angstrom':
            factor = 1
        else:
            raise ValueError('Unsupported unit: {}'.format(unit))

        if step >= 1:
            step -= 1

        idx_inporient = [i for i,s in enumerate(self.__logdata) if 'Input orientation:' in s][step]

        pedtab = Chem.GetPeriodicTable()
        symblist = []
        xyzlist = []
        numHyphenLine = 0
        for line in self.__logdata[idx_inporient:]:
            if re.fullmatch('^-+$', line):
                numHyphenLine += 1
            if numHyphenLine == 3:
                break
            if numHyphenLine > 2:
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


    def giveXYZBlock_in_InputOrientation(self, step=-1):
        symblist, xyzlist = self.giveMolecularGeometry_in_InputOrientation(step=step, unit='Angstrom')
        xyzblock = '\n'.join(['{} {} {} {}'.format(s,x,y,z) for s, (x,y,z) in zip(symblist,xyzlist)])
        return xyzblock

