import re as re
import sys as sys

#with open(sys.argv[1]) as f, open(sys.argv[2],'w') as output:
#	for line in f:
#		line=re.sub(r"(\:[0-9]+\.[0-9]+\:)",":",line)
#		output.write(line)	

with open(sys.argv[1]) as f, open(sys.argv[2],'w') as output:
	for line in f:
		line=re.sub(r"(.)*(\s)\({1}","(",line)
		output.write(line)
