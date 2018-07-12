#!/usr/bin/env python 
""" 
Read phonon dispersion from Gulp simulator
http://gulp.curtin.edu.au/gulp/
"""
#import .gulp_json
from Gulp_json_project.gulp_json import *
import argparse
import sys

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Read <prefix>.eig, <prefix>.gin and <prefix>.disp files from Gulp'
                                                 'and write a .json file to use in the phononwebsite.')
    parser.add_argument('prefix',           help='the prefix used in the .gin, .eig and .disp files calculation')
    parser.add_argument('-r','--reps',      help='number of default repetitions of the unit cell')
    parser.add_argument('-l','--labels',    help='string with the labels of the k-points. Eg. \"GMKG\" ')
    parser.add_argument('-w','--writeonly', help='only write json file (do not open the web browser)', action="store_true")

    #check if enough arguments are present
    if len(sys.argv)<2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    prefix = args.prefix
    
    q = gulp_json(prefix)
    if args.labels: q.set_labels(args.labels)
    if args.reps:   q.set_repetitions(args.reps)

    
    #q.write_json()
   
