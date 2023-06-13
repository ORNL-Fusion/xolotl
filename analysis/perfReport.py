#!/usr/bin/env python

import io
import math
import yaml

data = yaml.safe_load(io.open('/home/4pf/build/xolotl/Debug/perf.yaml', 'r'))
print(yaml.dump(data))

avgJacobianTime = data['Timers']['rhsJacobianTimer']['average']
print(avgJacobianTime)
