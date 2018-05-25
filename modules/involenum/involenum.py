import sys, re, pymongo, json
from itertools import combinations

client = pymongo.MongoClient("mongodb://manager:toric@129.10.135.170/ToricCY")
ToricCY = client.ToricCY

INVOL = ToricCY.INVOL

N = 6

divisors = {
    "dP": "{1,0,0,[1-9]}",
    "NSR": "{1,0,0,[1-9][^}]",
    "K3": "{1,0,1,20}",
    "Wilson": "{1,[1-9]+,0,[1-9]+}",
    "EWilson": "{1,1,0,[1-9]+}",
    "SD1": "{1,0,1,21}",
    "SD2": "{1,0,2,30}"
}

odims = {
    "O3": 3,
    "O5": 5,
    "O7": 7,
    "O9": 9
}

odims_keys = sorted(list(odims.keys()))

fields = [
    "INVOLS",
    "dP",
    "NSR",
    "K3",
    "Wilson",
    "EWilson",
    "SD1",
    "SD2",
    "dPK3",
    "dPWilson",
    "dPEWilson",
    "K3Wilson",
    "K3EWilson",
    "dPK3Wilson",
    "dPK3EWilson",
    "Free",
    "HO"
] + odims_keys \
  + [''.join(x) for x in list(combinations(odims_keys, 2))] \
  + [''.join(x) for x in list(combinations(odims_keys, 3))] \
  + [''.join(odims_keys)]

for line in iter(sys.stdin.readline, ''):
    geom_doc = json.loads(line.rstrip("\n"))

    invol_anomalies = {"O3O5": [], "O5O7": []}

    data_dict = {
        "TRIANGSR": {},
        "GEOMITENS": {},
        "GEOMITENSSR": {}
    }

    for field in fields:
        if field == "HO":
            data_dict["TRIANGSR"][field] = [0 for i in range(N - 1)]
            data_dict["GEOMITENS"][field] = [0 for i in range(N - 1)]
            data_dict["GEOMITENSSR"][field] = [0 for i in range(N - 1)]
        else:
            data_dict["TRIANGSR"][field] = 0
            data_dict["GEOMITENS"][field] = 0
            data_dict["GEOMITENSSR"][field] = 0

    geom_ntriangs_dict = {}
    geom_sr_ntriangs_dict = {}
    geom_data_dict = {}
    geom_sr_data_dict = {}

    geom_invol_anomalies = {"O3O5": [], "O5O7": [], "O3O9": [], "O7O9": []}

    has_invol_flag = False
    orientifold_done_flag = False

    invol_curs = INVOL.find({"H11": geom_doc['H11'], "POLYID": geom_doc['POLYID'], "GEOMN": geom_doc['GEOMN'], "OPLANES": {"$exists": True}}, \
                            {"_id": 0, "TRIANGN": 1, "INVOLN": 1, "INVOL": 1, "OPLANES": 1, "INVOLDIVCOHOM": 1, "ITENSXDINVOL": 1, "SRINVOL": 1, "SMOOTH": 1, "H11-": 1, "VOLFORMPARITY": 1}, \
                            no_cursor_timeout = True)
    for invol_doc in invol_curs:
        has_invol_flag = True
        if "OPLANES" in invol_doc:
            orientifold_done_flag = True
            nonzero_flag = False

            invol_dict = {}

            for field in fields:
                if field in divisors:
                    invol_dict[field] = sum(1 for cohom_string in invol_doc["INVOLDIVCOHOM"] if re.search(divisors[field], cohom_string))
                    if invol_dict[field] > 0:
                        nonzero_flag = True
                elif invol_doc["OPLANES"] and field in odims:
                    if (invol_doc["VOLFORMPARITY"] == -1 and field in ["O3", "O7"]) or (invol_doc["VOLFORMPARITY"] == 1 and field in ["O5", "O9"]):
                        invol_dict[field] = sum(1 for oplane in invol_doc["OPLANES"] if oplane["ODIM"] == odims[field])
                        if invol_dict[field] > 0:
                            nonzero_flag = True
                    else:
                        invol_dict[field] = 0
                elif field == "HO":
                    invol_dict[field] = [0 for i in range(N - 1)]
                else:
                    invol_dict[field] = 0

            if all(invol_dict[odim] == 0 for odim in odims):
                if invol_doc["SMOOTH"]:
                    invol_dict["Free"] = 1
                    nonzero_flag = True

            invol_dict["INVOLS"] = 1
            invol_dict["HO"][invol_doc["H11-"] - 1] = 1

            invol_dict["dPK3"] = invol_dict["dP"] * invol_dict["K3"]
            invol_dict["dPWilson"] = invol_dict["dP"] * invol_dict["Wilson"]
            invol_dict["dPEWilson"] = invol_dict["dP"] * invol_dict["EWilson"]
            invol_dict["K3Wilson"] = invol_dict["K3"] * invol_dict["Wilson"]
            invol_dict["K3EWilson"] = invol_dict["K3"] * invol_dict["EWilson"]
            invol_dict["dPK3Wilson"] = invol_dict["dP"] * invol_dict["K3"] * invol_dict["Wilson"]
            invol_dict["dPK3EWilson"] = invol_dict["dP"] * invol_dict["K3"] * invol_dict["EWilson"]

            invol_dict["O3O5"] = invol_dict["O3"] * invol_dict["O5"]
            invol_dict["O3O7"] = invol_dict["O3"] * invol_dict["O7"]
            invol_dict["O5O7"] = invol_dict["O5"] * invol_dict["O7"]
            invol_dict["O3O5O7"] = invol_dict["O3"] * invol_dict["O5"] * invol_dict["O7"]

            if invol_doc["ITENSXDINVOL"] and invol_doc["SRINVOL"]:
                for field in fields:
                    if field == "HO":
                        for i in range(N - 1):
                            data_dict["TRIANGSR"][field][i] += invol_dict[field][i]
                    else:
                        data_dict["TRIANGSR"][field] += invol_dict[field]

            if invol_doc["ITENSXDINVOL"]:
                if invol_doc['INVOL'] in geom_ntriangs_dict:
                    geom_ntriangs_dict[invol_doc['INVOL']] += 1                            
                else:
                    geom_ntriangs_dict[invol_doc['INVOL']] = 1
                    geom_data_dict[invol_doc['INVOL']] = {}
                    for field in fields:
                        if field == "HO":
                            geom_data_dict[invol_doc['INVOL']][field] = [0 for i in range(N - 1)]
                            for i in range(N - 1):
                                geom_data_dict[invol_doc['INVOL']][field][i] += invol_dict[field][i]
                        else:
                            geom_data_dict[invol_doc['INVOL']][field] = invol_dict[field]

                if invol_doc["SRINVOL"]:
                    if invol_doc['INVOL'] in geom_sr_ntriangs_dict:
                        geom_sr_ntriangs_dict[invol_doc['INVOL']] += 1
                    else:
                        geom_sr_ntriangs_dict[invol_doc['INVOL']] = 1
                        geom_sr_data_dict[invol_doc['INVOL']] = {}
                        for field in fields:
                            if field == "HO":
                                geom_sr_data_dict[invol_doc['INVOL']][field] = [0 for i in range(N - 1)]
                                for i in range(N - 1):
                                    geom_sr_data_dict[invol_doc['INVOL']][field][i] += invol_dict[field][i]
                            else:
                                geom_sr_data_dict[invol_doc['INVOL']][field] = invol_dict[field]

                if invol_dict["O3O5"]:
                    geom_invol_anomalies["O3O5"].append({"H11": geom_doc['H11'], "POLYID": geom_doc['POLYID'], "GEOMN": geom_doc['GEOMN'], "TRIANGN": invol_doc['TRIANGN'], "INVOLN": invol_doc['INVOLN'], "INVOL": invol_doc['INVOL'], "SRINVOL": invol_doc['SRINVOL']})

                if invol_dict["O5O7"]:
                    geom_invol_anomalies["O5O7"].append({"H11": geom_doc['H11'], "POLYID": geom_doc['POLYID'], "GEOMN": geom_doc['GEOMN'], "TRIANGN": invol_doc['TRIANGN'], "INVOLN": invol_doc['INVOLN'], "INVOL": invol_doc['INVOL'], "SRINVOL": invol_doc['SRINVOL']})

                if invol_dict["O3O9"]:
                    geom_invol_anomalies["O3O9"].append({"H11": geom_doc['H11'], "POLYID": geom_doc['POLYID'], "GEOMN": geom_doc['GEOMN'], "TRIANGN": invol_doc['TRIANGN'], "INVOLN": invol_doc['INVOLN'], "INVOL": invol_doc['INVOL'], "SRINVOL": invol_doc['SRINVOL']})

                if invol_dict["O7O9"]:
                    geom_invol_anomalies["O7O9"].append({"H11": geom_doc['H11'], "POLYID": geom_doc['POLYID'], "GEOMN": geom_doc['GEOMN'], "TRIANGN": invol_doc['TRIANGN'], "INVOLN": invol_doc['INVOLN'], "INVOL": invol_doc['INVOL'], "SRINVOL": invol_doc['SRINVOL']})

    for invol in geom_ntriangs_dict:
        if geom_ntriangs_dict[invol] == geom_doc['NTRIANGS']:
            for field in fields:
                if field == "HO":
                    for i in range(N - 1):
                        data_dict["GEOMITENS"][field][i] += geom_data_dict[invol][field][i]
                else:
                    data_dict["GEOMITENS"][field] += geom_data_dict[invol][field]

            for doc in geom_invol_anomalies["O3O5"]:
                if invol == doc["INVOL"]:
                    invol_anomalies["O3O5"].append(doc)

            for doc in geom_invol_anomalies["O5O7"]:
                if invol == doc["INVOL"]:
                    invol_anomalies["O5O7"].append(doc)

            for doc in geom_invol_anomalies["O3O9"]:
                if invol == doc["INVOL"]:
                    invol_anomalies["O3O9"].append(doc)

            for doc in geom_invol_anomalies["O7O9"]:
                if invol == doc["INVOL"]:
                    invol_anomalies["O7O9"].append(doc)

        if invol in geom_sr_ntriangs_dict and geom_sr_ntriangs_dict[invol] == geom_doc['NTRIANGS']:
            for field in fields:
                if field == "HO":
                    for i in range(N - 1):
                        data_dict["GEOMITENSSR"][field][i] += geom_sr_data_dict[invol][field][i]
                else:
                    data_dict["GEOMITENSSR"][field] += geom_sr_data_dict[invol][field]
    
    if has_invol_flag:
        data_dict.update({x: geom_doc[x] for x in ["H11", "POLYID", "GEOMN"]})
        data_dict.update(invol_anomalies)
        data_dict.update({"orientifoldMARK": orientifold_done_flag})

        print("set INVOLENUM " + json.dumps({x: geom_doc[x] for x in ["POLYID", "GEOMN"]}, separators = (',', ':')) + " " + json.dumps(data_dict, separators = (',', ':')))
    
    print("")
    sys.stdout.flush()

client.close()