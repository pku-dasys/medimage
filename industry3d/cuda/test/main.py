#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
data =''
iter=578

for i in (0.1, 0.01, 0.001, 0.0001):
    for j in (0.1, 0.01, 0.001, 0.0001):
        for k in (0.1, 0.01, 0.001, 0.0001):
            for m in (0.01, 0.001, 0.0001, 0.00001):
                for n in (0.01, 0.001, 0.0001, 0.00001):
					if i>0.001 or (i==0.001 and j==0.1):
						continue
					with open('test.json', 'r+') as f:
						for line in f.readlines():
							if (line.find('ALPHA') == 1):
								line = '"ALPHA" : ' + str(i) + ',\n'
							if (line.find('BETA') == 1):
								line = '"BETA" : ' + str(j) + ',\n'
							if (line.find('EPSILON') == 1):
								line = '"EPSILON" : ' + str(k) + ',\n'
							if (line.find('LAMBDA_IMG') == 1):
								line = '"LAMBDA_IMG" : ' + str(m) + ',\n'
							if (line.find('LAMBDA_EDGE') == 1):
								line = '"LAMBDA_EDGE" : ' + str(n) + ',\n'
							data += line
					with open('test.json', 'w') as f:
						f.writelines(data)
					os.system("../ct3d test.json")
					os.system("mv IMG IMG_" + str(iter))
					os.system("mv EDGE EDGE_" + str(iter))
					os.system("rm test.json")
					os.system("cp test_180x256x256_128_cone_shepp_logan_str.json test.json")
					iter=iter+1
					data = ''




