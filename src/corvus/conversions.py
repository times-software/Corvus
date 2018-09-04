isNum = lambda x: isinstance(x, (int, float)) 
isStr = lambda x: isinstance(x, basestring)
isList = lambda x: isinstance(x, list)

def bohr2angstrom(x):
    if isinstance(x, (int, float)):
        return x * 0.529177249
    elif isList(x) and all(map(isNum, x)):
        return [xi * 0.529177249 for xi in x]

def angstrom2bohr(x):
    if isinstance(x, (int, float)):
        return x / 0.529177249
    elif isList(x) and all(map(isNum, x)):
        return [xi / 0.529177249 for xi in x]
