import numpy as np
import sys as sys
import re

#f=file.open(sys.argv[1])
with open(sys.argv[1],'ra+') as f:
	with open(sys.argv[2],'a') as j:
		for line in f:
			line=re.sub(r"\[[0-9]*\]","",line)	
			j.write(line)

