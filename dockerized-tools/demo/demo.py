#!/usr/bin/env python

# Test script
# Simple example to demo docker


# INPUT:
# filenames      list of filenames

# OUTPUT:
# file_list.txt  text file containing list of filenames


# python script to align RNA-seq reads
VERSION = "0.0.1"

import argparse
    
def demo(args):
    
    f = open('file_list.txt', 'w')
    
    for filename in args.filenames:
        f.write(filename + '\n')
    
    f.close()
        
if __name__ == "__main__":
    """ Parse the command line arguments """
    parser = argparse.ArgumentParser(description='Report list of filenames')
    parser.add_argument('filenames', metavar='filenames', nargs='+', help='filenames including paths')
    parser.add_argument("--version", action='version', version=VERSION) 
    args = parser.parse_args()

    """ Run the desired methods """
    demo(args)




