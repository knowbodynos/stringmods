#!/shared/apps/sage-7.4/local/bin/sage -python

import sys, json
from time import sleep

for line in iter(sys.stdin.readline,''):
    testdoc=json.loads(line.rstrip("\n"))

    input = testdoc['INPUT']

    sleep(30)

    print("+TEST." + json.dumps({'INPUT': input}, separators = (',', ':')) + ">" + json.dumps({'OUTPUT': 0 - input}, separators = (',', ':')))
    print("@")
    sys.stdout.flush()