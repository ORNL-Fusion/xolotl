#!/usr/bin/env python

import glob
import io
import os
import sys
import math
import yaml
import matplotlib
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


def getDataFromFiles(rootDir):
    files = glob.glob(f"{baseName}perf_r[0-9]*.yaml", root_dir = rootDir)
    if not files:
        raise RuntimeError("No files found")

    data = []
    for f in files:
        print(f)
        d = yaml.safe_load(io.open(rootDir + '/' + f, 'r'))
        data.insert(d['Rank'], d)

    return data


def plotAll(ax1, ax2, data):
    makePlot(ax1, ['solveTimer', 'rhsFunctionTimer', 'Flux',
        'rhsJacobianTimer', 'Partial Derivatives'], 'Timers', data)
    makePlot(ax2, ['Flux', 'Partial Derivatives'], 'Counters', data)


def plotLocal(ax1, ax2, data):
    makePlot(ax1, ['rhsFunctionTimer', 'Flux',
        'rhsJacobianTimer', 'Partial Derivatives'], 'Timers', data)
    makePlot(ax2, ['Flux', 'Partial Derivatives'], 'Counters', data)


def getAverage(itemName, data):
    avg = 0.0
    for d in data:
        avg += d['Timers'][itemName]
    avg /= len(data1)
    return avg


baseName = ''
if len(sys.argv) > 1:
    baseName = sys.argv[1] + '_'

rootDir1 = os.getcwd()
if len(sys.argv) > 2:
    rootDir1 = sys.argv[2]

data1 = getDataFromFiles(rootDir1)

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

if len(sys.argv) > 3:
    fig, axes = plt.subplots(2, 2, gridspec_kw={'width_ratios': [3, 1]})
    fig.suptitle(baseName)
    rootDir2 = sys.argv[3]
    data2 = getDataFromFiles(rootDir2)
    # plotAll(axes[0][0], axes[0][1], data1)
    # plotAll(axes[1][0], axes[1][1], data2)
    plotLocal(axes[0][0], axes[0][1], data1)
    plotLocal(axes[1][0], axes[1][1], data2)
    print(getAverage('solveTimer', data1))
    print(getAverage('solveTimer', data2))
    print(getAverage('rhsJacobianTimer', data1))
    print(getAverage('rhsJacobianTimer', data2))
else:
    fig, axes = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})
    fig.suptitle(baseName)
    plotFromFiles(axes[0], axes[1], data1)
    print(getAverage('rhsJacobianTimer', data1))

plt.show()

