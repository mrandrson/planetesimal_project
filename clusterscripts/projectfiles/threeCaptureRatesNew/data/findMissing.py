import numpy as np
import os

def findMissing(directory):
	 a = os.listdir(directory)
	 for i in range(0, 15):
		 for j in range(0, 79):
			 s = 'v' + str(i) + 'b' + str(j) + '.txt'
			 if s not in a:
				 print(i, j)
				 
