
def __factor(unit):
    if unit in ['Angstrom', 'A', 'Ang']:
        return 0.529177210544
    elif unit in ['Bohr', 'a.u.']:
        return 1
    else:
        raise ValueError('invalid unit: {}'.format(unit))

def checkValidUnit(unit):
    return unit in ['Angstrom', 'A', 'Ang', 'Bohr', 'a.u.']

def getUnitConversionFactor(oldunit=None, newunit=None):
    if oldunit is None:
        raise ValueError('oldunit is None')
    if newunit is None:
        raise ValueError('newunit is None')

    oldfactor = __factor(oldunit)
    newfactor = __factor(newunit)

    return newfactor / oldfactor

