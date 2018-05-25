import sys, json

def mat2py(mat):
    py = mat.replace('{', '[').replace('}', ']')
    return eval(py)

def py2mat(py):
    mat = str(py).replace(' ', '').replace('[', '{').replace(']', '}')
    return mat

for line in iter(sys.stdin.readline, ''):
    triang_doc = json.loads(line.rstrip("\n"))

    query = {x: triang_doc[x] for x in ["POLYID", "GEOMN", "TRIANGN"]}
    output = {x: py2mat(triang_doc[x]) for x in ["TRIANG"]}

    print("set TRIANG " + json.dumps(query, separators = (',', ':')) + " " + json.dumps(output, separators = (',', ':')))
    print("")
    sys.stdout.flush()