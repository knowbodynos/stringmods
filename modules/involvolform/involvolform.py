import sys, re, json
from sympy.combinatorics import Permutation

def deldup(lst):
    "Delete duplicate elements in lst."
    return [lst[i] for i in range(len(lst)) if lst[i] not in lst[:i]]

def transpose_list(lst):
    "Get the transpose of a list of lists."
    try:
        if len(lst) == 0:
            raise IndexError('List is empty.')
        elif not all([len(x) == len(lst[0]) for x in lst]):
            raise IndexError('Lists have different sizes.')
        else:
            return [[lst[i][j] for i in range(len(lst))] for j in range(len(lst[0]))]
    except IndexError:
        return []

def mat2py(mat):
    py = mat.replace('{', '[').replace('}', ']')
    return eval(py)

def cycle(arr, orig = None):
    if arr == orig:
        return []
    if not orig:
        orig = arr
    new_arr = arr[:]
    new_arr.append(new_arr.pop(0))
    return [arr] + cycle(new_arr, orig = orig)

def is_same_parity(perm0, perm1):
    """Check if 2 permutations are of equal parity.

    Assume that both permutation lists are of equal length
    and have the same elements. No need to check for these
    conditions.
    """
    perm1 = perm1[:] ## copy this list so we don't mutate the original

    transCount = 0
    for loc in range(len(perm0) - 1):                         # Do (len - 1) transpositions
        p0 = perm0[loc]
        p1 = perm1[loc]
        if p0 != p1:
            sloc = perm1[loc:].index(p0) + loc          # Find position in perm1
            perm1[loc], perm1[sloc] = p0, p1          # Swap in perm1
            transCount += 1

    # Even number of transpositions means equal parity
    if transCount % 2 == 0:
        return True
    else:
        return False

def volform_parity(ndivs, invol):
    order = list(range(ndivs))
    for x, y in invol:
        temp = order[x]
        order[x] = order[y]
        order[y] = temp
    perm = Permutation(order)
    cycles = cycle(list(range(ndivs)))
    volform_term_parities = [is_same_parity(cycles[i][1:], [i^perm for i in cycles[i^perm][1:]]) for i in range(ndivs)]
    if all(x for x in volform_term_parities):
        return 1
    elif all(not x for x in volform_term_parities):
        return -1
    else:
        return 0

for line in iter(sys.stdin.readline,''):
    invol_doc = json.loads(line.rstrip("\n"))
    ndivs = invol_doc["H11"] + 4

    invol = deldup([sorted([y - 1 for y in x]) for x in mat2py(re.sub("D([0-9]+)->D([0-9]+)", r"[\1,\2]", invol_doc['INVOL']))])
    parity = volform_parity(ndivs, invol)

    invol_query = {key: invol_doc[key] for key in ["POLYID", "GEOMN", "TRIANGN", "INVOLN"]}
    print("set INVOL " + json.dumps(invol_query, separators=(',', ':')) + " " + json.dumps({"VOLFORMPARITY": parity}, separators = (',', ':')))
    print("")
    sys.stdout.flush()