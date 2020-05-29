#!/usr/bin/python

import glob
import os
import matplotlib.pyplot as plt

def main():
    for filename in glob.glob('*.xvg'):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            dist = []
            energy = []
            for line in f:
                dist.append(line.split()[0])
                energy.append(line.split()[1])
    
            dist.sort()
            plt.plot(dist, energy)
            
    plt.xlabel('distance')
    plt.ylabel('COM Pull energy')
    plt.show()
    

if __name__ == "__main__":
    main()
