#!/usr/bin/env python

import glob
import io
import math
import yaml
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

def makePlot(x, groupName, data):
    nranks = len(data)
    y = [ [] for _ in range(nranks) ]
    for i, d in enumerate(data):
        for name in x:
            y[i].append(d[groupName][name])

    nx = len(x)
    df = pd.DataFrame(
            zip(x*nranks, ["0"]*nx+["1"]*nx+["2"]*nx+["3"]*nx, y[0]+y[1]+y[2]+y[3]),
            columns=[groupName, "Rank", "Value"])
    plt.figure()
    sb.barplot(x=groupName, hue="Rank", y="Value", data=df)


rootDir='/home/4pf/build/xolotl/serial/'
files = glob.glob('perf_r[0-9]*.yaml', root_dir = rootDir)

data = []
for f in files:
    d = yaml.safe_load(io.open(rootDir + f, 'r'))
    data.insert(d['Rank'], d)

makePlot(['Flux', 'Partial Derivatives'], 'Timers', data)
makePlot(['Flux', 'Partial Derivatives'], 'Counters', data)

plt.show()

# for d in data:
#     print(yaml.dump(d))

