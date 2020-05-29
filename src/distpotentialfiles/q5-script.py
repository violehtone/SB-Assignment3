#!/usr/bin/python

import glob
import os
import matplotlib.pyplot as plt

def main():
    for filename in glob.glob('*.xvg'):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            distpotential = {}
            
            for line in f:
                distpotential[line.split()[0]] = line.split()[1]
            
            lists = sorted(distpotential.items())
            dist, potential = zip(*lists)
            plt.plot(dist, potential)
            
    plt.xlabel('distance')
    plt.ylabel('Potential energy')
    plt.show()
    

if __name__ == "__main__":
    main()
