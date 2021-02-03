
from scan_reader import MzXML
from extract import PeakExtractor, calculate_isolation_specificity

import csv
import re

def read_isolation_mz(filter_line):
    match = re.search(r'\s([\d\.]+)@', filter_line)
    if match:
        return float(match.group(1))
    return 0


def process_peptides(peptides_file_path, scans_file, quant_ini, output_file, verbose=True):
    isolation_window = 0.4

    extractor = PeakExtractor()
    extractor.load_options(quant_ini)

    if verbose:
        print("Read quant.ini")

    scans = MzXML()
    scans.open(scans_file)
    if verbose:
        print("opened MzXML file")

    with open(peptides_file_path) as peptides_file:
        with open(output_file, 'w', newline="") as outfile:
            peptides = csv.DictReader(peptides_file)

            quant_fields = []
            for field_type in ['Sn', 'Mz Error']:
                for channel in extractor.channel_names():
                    quant_fields.append(channel + ' ' + field_type)
            quant_fields += ['Ion Injection Time',
                             'Precursor Intensity',
                             'Elapsed Scan Time',
                             'Isolation Mz',
                             'Isolation Specificity',
                             'Sum Sn']

            fieldnames = peptides.fieldnames + quant_fields

            writer = csv.DictWriter(outfile, fieldnames)
            writer.writeheader()

            index = 0
            for peptide in peptides:
                # m scripts dont handle enclosures so scrub commas
                for key in peptide:
                    peptide[key] = peptide[key].replace(',', '.')

                scan_number = int(peptide['MS/MS scan number'])
                peptide_mz = float(peptide['Mass'])/float(peptide['Charge']) + 1.007276466812

                index += 1
                if verbose and index % 1000 == 0:
                    print("Working on peptide #" + str(index))
                
                scan = scans.get_scan(scan_number)
                
                extractor.extract_peaks(scan, peptide_mz)
                
                ion_injection_time = scan['ionInjectionTime']
                
                
                scan_time = scan['elapsedScanTime']
                precursor_intensity = scan['precursor_intensity']
                isolation_mz = read_isolation_mz(scan['filterLine'])

                parent_scan_number = scan['precursor_scan_num']
                parent_scan = scans.get_scan(parent_scan_number)
                isolation_specificity = calculate_isolation_specificity(parent_scan,
                                                                        isolation_mz,
                                                                        scan['precursor_charge'],
                                                                        isolation_window)

                sum_sn = 0
                for reporter in extractor.reporters:
                    key = reporter.name + ' Sn'
                    # TODO
                    noise = 1
                    if reporter.noise > 0:
                        noise = reporter.noise
                    sum_sn += reporter.intensity / noise
                    peptide[key] = reporter.intensity / noise

                for reporter in extractor.reporters:
                    key = reporter.name + ' Mz Error'
                    peptide[key] = reporter.mz_error

                peptide['Sum Sn'] = sum_sn
                peptide['Ion Injection Time'] = ion_injection_time
                peptide['Elapsed Scan Time'] = scan_time
                peptide['Precursor Intensity'] = precursor_intensity
                peptide['Isolation Mz'] = isolation_mz
                peptide['Isolation Specificity'] = isolation_specificity
                writer.writerow(peptide)