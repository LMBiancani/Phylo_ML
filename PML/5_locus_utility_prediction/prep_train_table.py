'''
A script to prepare a table for model training / prediction
Concatenates assessed properties tables from several datasets
and removes selected columns
'''
import sys
import csv
import os

def output_line(line, handle, columns_to_exclude):
	'''
	A function to exclude selected columns and NA vals
	'''
	if len(columns_to_exclude) > 0:
		for idx in columns_to_exclude:
			if idx < len(line):
				line.pop(idx)
	if "NA" not in line:
		rejoinedline = "\t".join(line)
		print(rejoinedline, file=handle)


def parse_path_file(pathfile, columns_to_exclude=[]):
	'''
	A function to parse input files and write output
	'''
	print ("read path file")
	working_dir=os.getcwd()
#	basename = os.path.basename(os.path.dirname(working_dir))
	with open(pathfile) as pathfilehandle:
		with open(sys.argv[2], "w") as outhandle:
			headerB = True
			for pathline in pathfilehandle:
				strippathline = pathline.strip()
				print(strippathline)
				with open(strippathline) as fhandle:
					reader = csv.reader(fhandle, delimiter='\t', quotechar='"')
					header = reader.__next__()
					if headerB:
						output_line(header, outhandle, columns_to_exclude)
						headerB = False
					for readerline in reader:
						readerline[0] = readerline[0]
						output_line(readerline, outhandle, columns_to_exclude)
	print("done")


def main():
	'''
	The main function
	'''
	if len(sys.argv) > 1:
		pathfile = sys.argv[1]
		if len(sys.argv) > 3:
			columns_to_exclude = sorted([int(x) for x in sys.argv[3:]], reverse=True)
			if 0 in columns_to_exclude:
				print ("should not be removing column 0")
			else:
				print ("excluding", columns_to_exclude)
				parse_path_file(pathfile, columns_to_exclude)
		else:
			parse_path_file(pathfile)
	else:
		print("python3 prep_train_table.py [file with a list of paths] [outfile] ([columns])")
		print("python3 prep_train_table.py train_list.txt outfile.txt")
		print("python3 prep_train_table.py train_list.txt outfile.txt 1 2 5")


if __name__ == "__main__":
	main()
