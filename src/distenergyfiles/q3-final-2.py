#!/usr/bin/python

import glob
import os
import matplotlib.pyplot as plt

def main():
    for filename in glob.glob('*.xvg'):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            distenergy = {}
            
            for line in f:
                distenergy[line.split()[0]] = line.split()[1]
            
            lists = sorted(distenergy.items())
            dist, energy = zip(*lists)
            plt.plot(dist, energy)
            
    plt.xlabel('distance')
    plt.ylabel('COM Pull energy')
    plt.show()
    

if __name__ == "__main__":
    main()
