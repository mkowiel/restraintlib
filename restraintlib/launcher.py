from __future__ import print_function

import argparse
import sys
import six

from restraintlib.printer import BusterPrinter
from restraintlib.printer import CsvPrinter
from restraintlib.printer import ShelxPrinter
from restraintlib.printer import PhenixPrinter
from restraintlib.printer import RefmacPrinter
from restraintlib.restraints import load_restraints_lib
from restraintlib.restraints import parse_pdb


class AllowedRestraintsConfig(object):

    def __init__(
        self,
        po4=True,
        po4terminal=True,
        nucleicacidbases=True,
        nucleicacidisobases=True,
        ribose_deoxyribose=True,
        ribose_deoxyribose_terminal=True,
    ):
        self.po4 = po4
        self.po4terminal = po4terminal
        self.nucleicacidbases = nucleicacidbases
        self.nucleicacidisobases = nucleicacidisobases
        self.ribose_deoxyribose = ribose_deoxyribose
        self.ribose_deoxyribose_terminal = ribose_deoxyribose_terminal


class RestraintLibLauncher(object):

    def __init__(self, log_stream=None):
        self.log_stream = sys.stdout
        if(type(log_stream) == str or type(log_stream) == six.text_type):
            self.log_stream = open(log_stream, 'w')
        elif log_stream is not None:
            self.log_stream = log_stream

    def report_progress(self, percent, total):
        pass

    def produce_restraints(self, in_pdb, printer, allowed_restraints_config):
        if printer == 'P':
            printer_cls = PhenixPrinter
        elif printer == 'S':
            printer_cls = ShelxPrinter
        elif printer == 'C':
            printer_cls = CsvPrinter
        elif printer == 'B':
            printer_cls = BusterPrinter
        else:
            printer_cls = RefmacPrinter

        restraint_list = load_restraints_lib()
        allowed_restraint_list = load_restraints_lib(
            allowed_restraints_config.po4,
            allowed_restraints_config.po4terminal,
            allowed_restraints_config.nucleicacidbases,
            allowed_restraints_config.nucleicacidisobases,
            allowed_restraints_config.ribose_deoxyribose,
            allowed_restraints_config.ribose_deoxyribose_terminal,
        )

        parse_pdb(in_pdb, restraint_list, allowed_restraint_list, self.log_stream, printer_cls)

        return []


def main():
    parser = argparse.ArgumentParser(description='Generate olgonucleotides restraints for pdb or mmcif file')
    parser.add_argument('printer', type=str, choices=['refmac', 'phenix', 'shelxl', 'buster', 'csv'], default='refmac',
                        help='Restraint output format')
    parser.add_argument('in_filename', type=str, default='in.pdb', help='Input file')
    parser.add_argument('out_filename', type=str, default='restraints.txt', help='Output restraints file')

    args = parser.parse_args()
    printer = args.printer.lower()
    in_pdb = args.in_filename
    out_filename = args.out_filename

    if printer == 'refmac':
        printer_cls = RefmacPrinter
    elif printer == 'phenix':
        printer_cls = PhenixPrinter
    elif printer == 'shelxl':
        printer_cls = ShelxPrinter
    elif printer == 'csv':
        printer_cls = CsvPrinter
    elif printer == 'buster':
        printer_cls = BusterPrinter
    else:
        print("Unknown printer {}, should be one of refmac, phenix, shelxl, buster, csv".format(printer))
        return

    restraint_list = load_restraints_lib()
    parse_pdb(in_pdb, restraint_list, restraint_list, out_filename, printer_cls)


if __name__ == "__main__":
    main()
