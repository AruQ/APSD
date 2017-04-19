#!/usr/bin/python

import matplotlib.pyplot as plt
import subprocess
import numpy as np
import sys


def myrun(cmd,nRows,nCol):
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		plt.axis([0, nRows-1, 0, nCol-1])
		plt.ion()
		A = []
		#stdout = []
		manager = plt.get_current_fig_manager()
		manager.resize(*manager.window.maxsize())
		step = 1
		while True:
				line = p.stdout.readline()
				#stdout.append(line)
				#print line,
				splitted = line.split( )
				#A = np.random.randint(2, size=(10, 10))
				splitted = map(int, splitted)
				A.append(splitted)
				if len(A) == nRows:
					plt.imshow(A)
					plt.pause(0.0001)
					A = []
					print step
					step += 1
				if line == '' and p.poll() != None:
						break
		#return ''.join(stdout)


if __name__ == "__main__":
	cmd = "./a.out " + sys.argv[1] +" "+sys.argv[2] +" "+ sys.argv[3]
	myrun(cmd, int(sys.argv[1]),int(sys.argv[2]))
