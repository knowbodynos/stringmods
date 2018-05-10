import sys, operator, itertools, json 

def mat2py(mat):
    py = mat.replace('{', '[').replace('}', ']')
    return eval(py)

def py2mat(py):
    mat = str(py).replace(' ', '').replace('[', '{').replace(']', '}')
    return mat

def disjointNIDs(NIDpairs, swaplist = []):
    if len(NIDpairs) == 0:
        return swaplist
    fullswaplist = []
    firstflag = True
    for NIDpair in NIDpairs:
        newNIDpairs = [x for x in NIDpairs if (x[0] > min(NIDpair)) and (x[1] > min(NIDpair)) and all([i not in x for i in NIDpair])]
        if len(swaplist) == 0:
            newswaplist = [[NIDpair]]
        else:
            newswaplist = [swaplist[-1] + [NIDpair]]
        if firstflag:
            newswaplist = swaplist + newswaplist
        fullswaplist += disjointNIDs(newNIDpairs, newswaplist)
        firstflag = False
    if len(swaplist) == 0:
        fullswaplist = sorted(fullswaplist, key = lambda x: (len(x), x[0]))
    return fullswaplist

for line in iter(sys.stdin.readline,''):
    triang_doc = json.loads(line.rstrip("\n"))

    rescws = mat2py(triang_doc['RESCWS'])
    divcohom = mat2py(triang_doc['DIVCOHOM'])
    itensXD = mat2py(triang_doc['ITENSXD'])
    SRideal = triang_doc['SRIDEAL']
    SRsets = sorted([[y - 1 for y in eval(("[" + x + "]").replace("D", "").replace("*", ","))] for x in SRideal.lstrip("{").rstrip("}").split(",")], key = lambda x: (len(x), operator.itemgetter(*range(len(x)))(x)))
    itensXDsets = []
    for i in range(len(rescws)):
        for j in range(i, len(rescws)):
            for k in range(j, len(rescws)):
                itensXDsets.append([[i, j, k], itensXD[i][j][k]])
    itensXDsets = sorted(itensXDsets, key = lambda x: (len(x[0]), operator.itemgetter(*range(len(x[0])))(x[0])))

    newdivcohom = [py2mat(x) for x in divcohom]

    NIDpairs = []
    for i in range(len(rescws)):
        for j in range(i + 1, len(rescws)):
            if (divcohom[i] == divcohom[j]) and (rescws[i] != rescws[j]):
                NIDpairs.append([i, j])

    #disjointsets = []
    #for i in range(len(NIDpairs)):
    #    combs = list(itertools.combinations(NIDpairs, i + 1))
    #    for comb in combs:
    #        if all([not any([k in y for k in comb[j] for y in comb[:j] + comb[j + 1:]]) for j in range(len(comb))]):
    #            disjointsets.append(list(comb))

    disjointsets = disjointNIDs(NIDpairs)

    allowedinvols = []
    involn = 1
    for invol in disjointsets:
        newSRsets = []
        for SRset in SRsets:
            newSRset = SRset
            for x in invol:
                newSRset = [x[1] if y == x[0] else x[0] if y == x[1] else y for y in newSRset]
            newSRset = sorted(newSRset)
            newSRsets.append(newSRset)
        newSRsets = sorted(newSRsets, key = lambda x: (len(x), operator.itemgetter(*range(len(x)))(x)))
        
        newitensXDsets = []
        for itensXDset in itensXDsets:
            newitensXDset = itensXDset[0]
            for x in invol:
                newitensXDset = [x[1] if y == x[0] else x[0] if y == x[1] else y for y in newitensXDset]
            newitensXDset = [sorted(newitensXDset), itensXDset[1]]
            newitensXDsets.append(newitensXDset)
        newitensXDsets = sorted(newitensXDsets, key = lambda x: (len(x[0]), operator.itemgetter(*range(len(x[0])))(x[0])))
        
        if (newSRsets == SRsets) or (newitensXDsets == itensXDsets):
            matinvol = "{" + ",".join([",".join(["D" + str(x[0] + 1) + "->D" + str(x[1] + 1), "D" + str(x[1] + 1) + "->D" + str(x[0] + 1)]) for x in invol]) + "}"
            SRinvol = (newSRsets == SRsets)
            itensXDinvol = (newitensXDsets == itensXDsets)
            involquery = {"POLYID": triang_doc['POLYID'], "GEOMN": triang_doc['GEOMN'], "TRIANGN": triang_doc['TRIANGN'], "INVOLN": involn}
            newinvolout = {"H11": triang_doc['H11'], "POLYID": triang_doc['POLYID'], "GEOMN": triang_doc['GEOMN'], "TRIANGN": triang_doc['TRIANGN'], "INVOLN": involn, "INVOL": matinvol, "INVOLDIVCOHOM": [py2mat(divcohom[x[0]]) for x in invol], "SRINVOL": SRinvol, "ITENSXDINVOL": itensXDinvol}
            print("set INVOL " + json.dumps(involquery, separators = (',', ':')) + " " + json.dumps(newinvolout, separators = (',', ':')))
            sys.stdout.flush()
            involn += 1
    involn -= 1

    triangquery = {"POLYID": triang_doc['POLYID'], "GEOMN": triang_doc['GEOMN'], "TRIANGN": triang_doc['TRIANGN']}
    newtriangoutput = {"DIVCOHOM1": newdivcohom, "NINVOL": involn}

    print("set TRIANG " + json.dumps(triangquery, separators = (',', ':')) + " " + json.dumps(newtriangoutput, separators = (',', ':')))
    print("")
    sys.stdout.flush()