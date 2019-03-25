import argparse
import sys

parser = argparse.ArgumentParser(prog='bestcase-phasing.py', description=__doc__)
parser.add_argument('readlen', metavar='READLENGTH', help='theoretical read length')
args = parser.parse_args()

block_id = -1
previous_chromosome = None
previous_position = 0

for line in sys.stdin:
	if line[0] == '#':
		print(line[:-1])
		continue
	fields = line.split('\t')
	chrom = fields[0]
	position = int(fields[1])
	format_field = fields[8].rstrip().split(':')
	sample_field = fields[9].rstrip().split(':')
	gt_pos = format_field.index('GT')
	gt = sample_field[gt_pos]
	if gt == '.':
		continue
	if gt[0] == gt[2]:
		print(line[:-1])
		continue
	if previous_chromosome != chrom or (position - previous_position) > int(args.readlen):
		block_id += 1
	# assume PS tag is at the end ...
	if 'PS' in format_field:
		sample_field[-1] = str(block_id)
	else:
		format_field.append('PS')
		sample_field.append(str(block_id))
	sample_field[gt_pos] = '|'.join([gt[0],gt[2]])
	fields[8] = ':'.join(format_field)
	fields[9] = ':'.join(sample_field)
	print('\t'.join(fields))
	previous_chromosome = chrom
	previous_position = position
