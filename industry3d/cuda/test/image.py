#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from PIL import Image
import numpy as np

img=[]

with open('IMG_475','r+') as f:
	for line in f.readlines():
		line=line.split()
		row=[]
		for i in line:
			row.append(float(i))
		img.append(row)

img=np.array(img)
image=Image.fromarray(img)
image.show()
