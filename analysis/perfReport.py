#!/usr/bin/env python

import glob
import io
import os
import sys
import math
import yaml
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

def makePlot(axis, x, groupName, data):
    nranks = len(data)
    y = [ [] for _ in range(nranks) ]
    for i, d in enumerate(data):
        for name in x:
            y[i].append(d[groupName][name])

    nx = len(x)
    rankArray = []
    valArray = []
    for r in range(nranks):
        rankArray += [str(r)]*nx
        valArray += y[r]
    df = pd.DataFrame(zip(x*nranks, rankArray, valArray),
            columns=[groupName, "Rank", "Value"])
    sb.barplot(x=groupName, hue="Rank", y="Value", data=df, ax=axis)

baseName = ''
if len(sys.argv) > 1:
    baseName = sys.argv[1] + '_'

rootDir = os.getcwd()
if len(sys.argv) > 2:
    rootDir = sys.argv[2]

files = glob.glob(f"{baseName}perf_r[0-9]*.yaml", root_dir = rootDir)
if not files:
    raise RuntimeError("No files found")

data = []
for f in files:
    d = yaml.safe_load(io.open(rootDir + '/' + f, 'r'))
    data.insert(d['Rank'], d)

fig, axes = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})
fig.suptitle(baseName)

makePlot(axes[0], ['solveTimer', 'rhsFunctionTimer', 'Flux',
    'rhsJacobianTimer', 'Partial Derivatives'], 'Timers', data)
makePlot(axes[1], ['Flux', 'Partial Derivatives'], 'Counters', data)

plt.show()

# for d in data:
#     print(yaml.dump(d))

