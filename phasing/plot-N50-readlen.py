import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(prog='plot-N50-readlen.py', description=__doc__)
parser.add_argument('tsv', metavar='FILE', nargs='+', help='tsv files containing N50s')
parser.add_argument('output', metavar='OUTPUT', help='name of output file')
args = parser.parse_args()

def extract_readlen(filename):
	splitted = filename.split('.')
	for element in splitted:
		if element[0:9] == 'bestcase-':
			return int(element[9:])
	assert(False)

data = []

for filename in args.tsv:
	value = -1
	readlen = extract_readlen(filename)
	for line in open(filename, 'r'):
		splitted = line.split()
		if splitted[1] == 'ALL':
			value = int(splitted[21])
			continue
	assert(value is not -1)
	data.append( (readlen, value) )

# sort values low coverage -> high coverage
data = sorted(data, key=lambda x: x[0])
print(data)
readlengths = [i[0] for i in data]
n50s = [i[1] for i in data]

print('readlengths:', readlengths)
print('n50s:', n50s)

plt.plot(readlengths, n50s, 'bo-')
plt.xlabel('readlength [bp]')
plt.ylabel('N50 [bp]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig(args.output)
