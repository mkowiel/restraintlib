from __future__ import print_function

import sys
import six

from .printer import CsvPrinter
from .printer import ShelxPrinter
from .printer import PhenixPrinter
from .printer import RefmacPrinter

from .restraints import load_restraints_lib
from .restraints import parse_pdb


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
    if len(sys.argv) > 3:
        printer = sys.argv[1]
        in_pdb = sys.argv[2]
        out_filename = sys.argv[3]
    elif len(sys.argv) > 2:
        printer = 'Refmac'
        in_pdb = sys.argv[1]
        out_filename = sys.argv[2]
    else:
        print('usage: restraints.py [printer] in.pdb restraints.txt')
        print('       Printer is one of Refmac, Phenix, Shelxl, Csv. Default=Refmac')
        printer = 'Refmac'
        in_pdb = 'in.pdb'
        out_filename = 'restraints.txt'
    printer = printer.lower()

    if printer == 'refmac':
        printer_cls = RefmacPrinter
    elif printer == 'phenix':
        printer_cls = PhenixPrinter
    elif printer == 'shelxl':
        printer_cls = ShelxPrinter
    elif printer == 'csv':
        printer_cls = CsvPrinter
    else:
        print("Unknown printer {}, should be one of Refmac, Phenix, Shelxl".format(printer))
        return

    restraint_list = load_restraints_lib()
    parse_pdb(in_pdb, restraint_list, restraint_list, out_filename, printer_cls)



