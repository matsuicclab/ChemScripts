import re


class Input():
    """
    loader for GAMESS input
    """
    def __init__(self, **args):
        if 'filePath' in args.keys():
            self.__init__fromFile(**args)
        elif 'logData' in args.keys():
            self.__init__fromInpData(**args)
        else:
            raise ValueError('args must contain filePath or cubeData')

    def __init__fromFile(self, filePath=None):
        """
        load log file
        """
        # ファイル読み込み
        with open(filePath, mode='r') as f:
            inpData = [s.strip('\n') for s in f.readlines()]
        self.__init__fromInpData(inpData=inpData)

    def __init__fromInpData(self, inpData=None):
        def decompose(inpData):
            # キーの位置を特定
            inpData = re.sub(r'( \$[a-zA-Z0-9]+( |\n))', r' \1 ', inpData)
            matchResults = []
            for i,m in enumerate(re.finditer(r' \$[a-zA-Z0-9]+( |\n)', inpData)):
                startIdx = m.start()
                endIdx = m.end()
                key = m.group().strip()
                if key == '$END' and i % 2 == 0:
                    raise RuntimeError()
                matchResults.append([i,startIdx,endIdx,key])
            # 各キーと値を取得
            d = {}
            for (_,start1,end1,key1),(_,start2,end2,key2) in zip(matchResults[0::2], matchResults[1::2]):
                key = re.sub(r'\$', '', key1)
                value = inpData[end1:start2]
                # コメント文を除去
                value_trimmedcomment = []
                for line in value.split('\n'):
                    for i, ch in enumerate(line):
                        if ch == '!':
                            count_s = line[:i].count('"')
                            count_d = line[:i].count("'")
                            if count_s % 2 == 0 and count_d % 2 == 0:
                                line = line[:i]
                                break
                    value_trimmedcomment.append(line)
                value = '\n'.join(value_trimmedcomment)
                if '=' in value:
                    d[key] = dict()
                    decomposedvalue = re.sub(r'=', r'\n=\n', re.sub(r'([^, ]+=)', r'\n\1', value)).split()
                    equalIndices = [i for i, line in enumerate(decomposedvalue) if line=='=']
                    for i,idx in enumerate(equalIndices):
                        childKeyIdx = idx-1
                        nextChildKeyIdx = len(decomposedvalue) if (i+1==len(equalIndices)) else equalIndices[i+1]-1
                        childKey = decomposedvalue[childKeyIdx]
                        childValueIdx = idx+1
                        childValue = re.sub(r'[, ]+$', '', ','.join(decomposedvalue[childValueIdx:nextChildKeyIdx]))
                        d[key][childKey] = childValue
                else:
                    d[key] = value
            return d

        def convertType(d):
            def __convert(v):
                try:
                    _v = int(v)
                except:
                    try:
                        _v = float(re.sub('d','e', v, 1))
                    except:
                        if v.startswith('"') and v.endswith('"'):
                            _v = v[1:-1]
                        elif v.startswith("'") and v.endswith("'"):
                            _v = v[1:-1]
                        else:
                            _v = v
                return _v

            for key, value in d.items():
                if type(value) is dict:
                    convertType(value)
                elif value == '.TRUE.':
                    d[key] = True
                elif value == '.FALSE.':
                    d[key] = False
                elif re.search(r'\([0-9]+\)', key):
                    d[key] = [__convert(v) for v in re.split(r'[, ]+', value)]
                else:
                    d[key] = __convert(value)
            return d
        #-------------------------------------------

        self.__inpData = inpData
        inpDict = decompose('\n'.join(inpData))
        self.__inpDict = convertType(inpDict)





