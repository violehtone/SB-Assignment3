#!/usr/bin/python

import argparse

def main(dist_file, energy_file):
    #create a new .xvg file with distance and potential combined for specific d
    with open('distpotential.xvg', 'w') as distpotentialfile:
        with open(dist_file, 'r') as distfile:
            with open(energy_file, 'r') as energyfile:
                for lined, linee in zip(distfile, energyfile):
                    if(lined[0] is not '#' and lined[0] is not '@' and linee[0] is not '#' and linee[0] is not '@'):
                        print(lined.split()[1], linee.split()[2], file=distpotentialfile)
                    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dist_file")
    parser.add_argument("-e", "--energy_file")
    args = parser.parse_args()

    dist_file = args.dist_file
    energy_file = args.energy_file
    main(dist_file, energy_file)
