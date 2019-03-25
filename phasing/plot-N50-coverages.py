import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(prog='plot-N50-coverages.py', description=__doc__)
parser.add_argument('tsv', metavar='FILE', nargs='+', help='tsv files containing N50s')
parser.add_argument('output', metavar='OUTPUT', help='name of output file')
args = parser.parse_args()

def extract_rate(filename):
	splitted = filename.split('.')
	for element in splitted:
		if element[0:4] == 'rate':
			return int(element[4:])
	assert(False)

data = []

for filename in args.tsv:
	value = -1
	rate = extract_rate(filename)
	for line in open(filename, 'r'):
		splitted = line.split()
		if splitted[1] == 'ALL':
			value = int(splitted[21])
			continue
	assert(value is not -1)
	data.append( (rate, value) )

# sort values low coverage -> high coverage
data = sorted(data, key=lambda x: x[0])
print(data)
coverages = [0.01*i[0] for i in data]
n50s = [i[1] for i in data]

print('coverages:', coverages)
print('n50s:', n50s)

plt.plot(coverages, n50s, 'bo-')
plt.xlabel('rate')
plt.ylabel('N50')
plt.savefig(args.output)
