#!/usr/bin/python

import os
import argparse
import matplotlib.pyplot as plt


def main(distenergy_file):
    dist = []
    energy = []
    with open(distenergy_file, 'r') as distenergyfile:
        for line in distenergyfile:
            dist.append(line.split()[0])
            energy.append(line.split()[1])
    
    dist.sort()
    plt.plot(dist, energy)
    plt.xlabel('distance')
    plt.ylabel('COM Pull energy')
    plt.show()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--distenergy_file")
    args = parser.parse_args()
    distenergy_file = args.distenergy_file
    main(distenergy_file)
