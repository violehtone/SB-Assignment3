#!/usr/bin/python

import argparse

def main(dist_file, energy_file):
    #create a new .xvg file with distance and energy combined for specific d
    with open('distenergy.xvg', 'w') as distenergyfile:
        with open(dist_file, 'r') as distfile:
            with open(energy_file, 'r') as energyfile:
                for lined, linee in zip(distfile, energyfile):
                    if(lined[0] is not '#' and lined[0] is not '@' and linee[0] is not '#' and linee[0] is not '@'):
                        print(lined.split()[1], linee.split()[1], file=distenergyfile)
                    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dist_file")
    parser.add_argument("-e", "--energy_file")
    args = parser.parse_args()

    dist_file = args.dist_file
    energy_file = args.energy_file
    main(dist_file, energy_file)
