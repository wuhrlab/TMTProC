import csv
import reporter_quant
import os
from sys import exit, argv
import argparse

iso_window_file = 'example_data/iso_window.txt'
do_quant = True

parser = argparse.ArgumentParser(description='TMTc Command line')
parser.add_argument('--peptides', required=True, help='The path to the csv formatted file with the list of peptides.')
parser.add_argument('--scans', required=True, help='The path to the mzXML file.')
parser.add_argument('--quant', required=True, help='An ini formatted file with the quant configuration.')
args = parser.parse_args(argv[1:])

peptides_file_path = args.peptides
scans_file = args.scans
quant_ini = args.quant

if not os.path.isfile(args.peptides):
    print("Peptides file does not exist: " + peptides_file_path)
    exit(1)
if not os.path.isfile(scans_file):
    print("Scans file does not exist: " + scans_file)
    exit(1)
if not os.path.isfile(quant_ini):
    print("Ini file does not exist: " + quant_ini)
    exit(1)

print("start")

pieces = peptides_file_path.split('.')
intermediate_file = pieces[0] + '_intermediate.csv'

if do_quant:
    print("Extracting peaks for peptides")
    reporter_quant.process_peptides(peptides_file_path, scans_file, quant_ini, intermediate_file)
    print("Done extracting peaks")

print("Initial quant done")
