#!/usr/bin/python

import os
import argparse
import matplotlib.pyplot as plt


def main(distenergy_file):
    distenergy = {}
    with open(distenergy_file, 'r') as distenergyfile:
        for line in distenergyfile:
            distenergy[line.split()[0]] = line.split()[1]

    lists = sorted(distenergy.items())
    dist, energy = zip(*lists)
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
