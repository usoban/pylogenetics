import csv
from numpy import zeros

def import_mega_csv(csvPath):
	"""
	Imports MEGA's CSV distance matrix and OTU names
	
	Returns a tuple (name[], distMtx[][])
	"""
	csvFile = open(csvPath)
	data    = csv.reader(csvFile, delimiter = ',')
	names   = []
	lines   = []
		
	for row in data:
		if len(row) > 0:
			names.append(row[0])
			lines.append(row[1:])
		else:
			break

	# First col is OTU name
	nspecies = len(lines[-1])
	distMtx  = zeros((nspecies, nspecies))
	distMtx[distMtx==0] = float('inf')
		
	for i in range(0, len(lines)):	
		for j in range(0, len(lines[i])):
			try:
				distMtx[i][j] = float(lines[i][j])
			except ValueError:
				break
				
	return (names, distMtx)
