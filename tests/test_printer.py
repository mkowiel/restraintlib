# -*- coding: utf-8 -*-
from unittest import TestCase

from restraintlib.atom import Atom
from restraintlib.printer import CsvPrinter
from restraintlib.printer import BusterPrinter
from restraintlib.printer import PhenixPrinter
from restraintlib.printer import PrinterBase
from restraintlib.printer import RefmacPrinter
from restraintlib.printer import ShelxPrinter
from restraintlib.printer import TuplePrinter
from restraintlib.restraints import Restraint
from restraintlib.restraints import load_restraints_lib


class PrinterTestCase(TestCase):

    def setUp(self):
        atoms = [
            Atom('A', '100', 'DT', 'OP1', 'A', (1.0, 0.0, 0.0), 1),
            Atom('A', '100', 'DT', 'P', 'A', (1.0, 1.6, 0.0), 2),
            Atom('A', '100', 'DT', 'OP2', 'A', (1.0, 3.2, 0.0), 3),
        ]

        self.restraints = [
            Restraint("angle", atoms, 117.6, 1.2, None, None, "aO1O2", "PO4==AA_0"),
            # for angle we need both lengths
            Restraint("dist", atoms[0:2], 1.487, 0.01, None, None, "dO1P4", "PO4==AA_0"),
            Restraint("dist", atoms[1:3], 1.483, 0.01, None, None, "dO2P4", "PO4==AA_0"),
        ]

        self.allowed_restraints = load_restraints_lib(
            po4=True,
            po4terminal=False,
            bases=False,
            isobases=False,
            ribose_deoxyribose=False,
            ribose_deoxyribose_terminal=False,
        )

        self.printer = PrinterBase(override_sigma=False)

        self.expected_lines_angle = [
            '',
        ]

        self.expected_lines_dist = [
            '',
        ]

    def test_get_angle(self):
        angle_text = self.printer.get_angle(self.restraints[0], self.restraints)

        for line in self.expected_lines_angle:
            self.assertIn(line, angle_text)

    def test_get_dist(self):
        dist = self.printer.get_dist(self.restraints[1], self.restraints)

        for line in self.expected_lines_dist:
            self.assertIn(line, dist)

    def test_print_restraints(self):
        restraints_lines = self.printer.print_restraints(self.restraints, self.allowed_restraints)

        expected_lines = self.expected_lines_angle + self.expected_lines_dist

        for line in expected_lines:
            self.assertIn(line, restraints_lines)

    def test_header(self):
        self.assertEqual(self.printer.header(), "")

    def test_footer(self):
        self.assertEqual(self.printer.footer(), "")


class ShelxPrinterTestCase(PrinterTestCase):

    def setUp(self):
        super(ShelxPrinterTestCase, self).setUp()
        self.printer = ShelxPrinter(override_sigma=False)

        self.expected_lines_angle = [
            'REM Chain A and resid 100 and name OP1 and altloc A renamed to OP1_A:100^a',
            'REM Chain A and resid 100 and name OP2 and altloc A renamed to OP2_A:100^a',
            'REM Restraint PO4==AA_0',
            'REM angle aO1O2 117.6 1.2',
            'DANG 2.540 0.016 OP1_A:100^a OP2_A:100^a',
        ]

        self.expected_lines_dist = [
            'REM Chain A and resid 100 and name OP1 and altloc A renamed to OP1_A:100^a',
            'REM Chain A and resid 100 and name P and altloc A renamed to P_A:100^a',
            'REM Restraint PO4==AA_0',
            'REM dist dO1P4 1.487 0.010',
            'DFIX 1.487 0.010 OP1_A:100^a P_A:100^a',
        ]


class PhenixPrinterTestCase(PrinterTestCase):

    def setUp(self):
        super(PhenixPrinterTestCase, self).setUp()
        self.printer = PhenixPrinter(override_sigma=False)

        self.expected_lines_angle = [
            '# Restraint PO4==AA_0 angle aO1O2 117.6 1.2',
            'angle {',
            'action = *change',
            'atom_selection_1 = chain A and resid 100 and name OP1 and altloc A',
            'atom_selection_2 = chain A and resid 100 and name P and altloc A',
            'atom_selection_3 = chain A and resid 100 and name OP2 and altloc A',
            'angle_ideal = 117.6',
            'sigma = 1.2',
            '}',
        ]

        self.expected_lines_dist = [
            '# Restraint PO4==AA_0 dist dO1P4 1.487 0.010',
            'bond {',
            'action = *change',
            'atom_selection_1 = chain A and resid 100 and name P and altloc A',
            'atom_selection_2 = chain A and resid 100 and name OP1 and altloc A',
            'distance_ideal = 1.487',
            'sigma = 0.010',
            'slack = None',
            '}',
        ]

    def test_header(self):
        self.assertEqual(self.printer.header(), "refinement.geometry_restraints.edits {")

    def test_footer(self):
        self.assertEqual(self.printer.footer(), "}")


class PhenixPrinterSigmaOverrideTestCase(PrinterTestCase):

    def setUp(self):
        super(PhenixPrinterSigmaOverrideTestCase, self).setUp()
        self.printer = PhenixPrinter(override_sigma=True)

        self.expected_lines_angle = [
            '# Restraint PO4==AA_0 angle aO1O2 117.6 1.2',
            'angle {',
            'action = *change',
            'atom_selection_1 = chain A and resid 100 and name OP1 and altloc A',
            'atom_selection_2 = chain A and resid 100 and name P and altloc A',
            'atom_selection_3 = chain A and resid 100 and name OP2 and altloc A',
            'angle_ideal = 117.6',
            'sigma = 3.0',
            '}',
        ]

        self.expected_lines_dist = [
            '# Restraint PO4==AA_0 dist dO1P4 1.487 0.010',
            'bond {',
            'action = *change',
            'atom_selection_1 = chain A and resid 100 and name P and altloc A',
            'atom_selection_2 = chain A and resid 100 and name OP1 and altloc A',
            'distance_ideal = 1.487',
            'sigma = 0.020',
            'slack = None',
            '}',
        ]

    def test_header(self):
        self.assertEqual(self.printer.header(), "refinement.geometry_restraints.edits {")

    def test_footer(self):
        self.assertEqual(self.printer.footer(), "}")


class PhenixPrinterO3C3TestCase(PrinterTestCase):

    def setUp(self):
        super(PhenixPrinterO3C3TestCase, self).setUp()
        self.printer = PhenixPrinter(override_sigma=True)

        atoms = [
            Atom('A', '100', 'T', "O3'", 'A', (1.0, 0.0, 0.0), 1),
            Atom('A', '100', 'T', "C3'", 'A', (1.0, 1.6, 0.0), 2),
            Atom('A', '100', 'T', "C4'", 'A', (1.0, 3.2, 0.0), 3),
        ]

        self.allowed_restraints = load_restraints_lib(
            po4=True,
            po4terminal=False,
            bases=False,
            isobases=False,
            ribose_deoxyribose=True,
            ribose_deoxyribose_terminal=False,
        )

        self.restraints = [
            Restraint("angle", atoms, 110.3, 2.4, None, None, "aO3'C3'C4'", "ribose_pyrimidine==Chi=anti"),
            # for angle we need both lengths
            Restraint("dist", atoms[0:2], 1.419, 0.009, None, None, "dO3'C3'", "ribose_pyrimidine==Chi=anti"),
            Restraint("dist", atoms[1:3], 1.525, 0.008, None, None, "dC3'C4'", "ribose_pyrimidine==Chi=anti"),
        ]

        self.expected_lines_angle = [
            "# Restraint ribose_pyrimidine==Chi=anti angle aO3'C3'C4' 110.3 2.4",
            'angle {',
            'action = *change',
            "atom_selection_1 = chain A and resid 100 and name O3' and altloc A",
            "atom_selection_2 = chain A and resid 100 and name C3' and altloc A",
            "atom_selection_3 = chain A and resid 100 and name C4' and altloc A",
            'angle_ideal = 110.3',
            'sigma = 3.0',
            '}',
        ]

        self.expected_lines_dist = [
            "# Restraint ribose_pyrimidine==Chi=anti dist dO3'C3' 1.419 0.009",
            'bond {',
            'action = *change',
            "atom_selection_1 = chain A and resid 100 and name O3' and altloc A",
            "atom_selection_2 = chain A and resid 100 and name C3' and altloc A",
            'distance_ideal = 1.419',
            'sigma = 0.020',
            'slack = None',
            '}',
        ]

    def test_header(self):
        self.assertEqual(self.printer.header(), "refinement.geometry_restraints.edits {")

    def test_footer(self):
        self.assertEqual(self.printer.footer(), "}")


class RefmacPrinterTestCase(PrinterTestCase):

    def setUp(self):
        super(RefmacPrinterTestCase, self).setUp()
        self.printer = RefmacPrinter(override_sigma=False)

        self.expected_lines_angle = [
            '# Restraint PO4==AA_0 angle aO1O2 117.6 1.2',
            'exte angle first chain A resi 100 atom OP1 altecode A next chain A resi 100 atom P altecode A',
            'next chain A resi 100 atom OP2 altecode A value 117.6 sigma 1.2 type 0',
        ]

        self.expected_lines_dist = [
            '# Restraint PO4==AA_0 dist dO1P4 1.487 0.010',
            'exte dist first chain A resi 100 atom OP1 altecode A second chain A resi 100 atom P altecode A',
            'value 1.487 sigma 0.010 type 0',
        ]


class CsvPrinterTestCase(PrinterTestCase):

    def setUp(self):
        super(CsvPrinterTestCase, self).setUp()
        self.printer = CsvPrinter(override_sigma=False)

        self.expected_lines_angle = [
            'angle,PO4==AA_0,aO1O2,A,100,OP1,A,A,100,P,A,A,100,OP2,A,117.6,1.2',
        ]

        self.expected_lines_dist = [
            'dist,PO4==AA_0,dO1P4,A,100,OP1,A,A,100,P,A,,,,,1.487,0.010',
        ]

    def test_header(self):
        self.assertNotEqual(self.printer.header(), "")

    def test_footer(self):
        self.assertEqual(self.printer.footer(), "")


class BusterPrinterTestCase(PrinterTestCase):

    def setUp(self):
        super(BusterPrinterTestCase, self).setUp()
        self.printer = BusterPrinter(override_sigma=False)

        self.expected_lines_angle = [
            '# Restraint PO4==AA_0 angle aO1O2 117.6 1.2',
            'NOTE BUSTER_UTILANGLE 117.6 1.2 A|100:OP1.A A|100:P.A A|100:OP2.A'
        ]

        self.expected_lines_dist = [
            '# Restraint PO4==AA_0 dist dO1P4 1.487 0.010',
            'NOTE BUSTER_DISTANCE =1.487 0.010 A|100:OP1.A A|100:P.A'
        ]


class BusterPrinterSigmaOverrideTestCase(PrinterTestCase):

    def setUp(self):
        super(BusterPrinterSigmaOverrideTestCase, self).setUp()
        self.printer = BusterPrinter(override_sigma=True)

        self.expected_lines_angle = [
            '# Restraint PO4==AA_0 angle aO1O2 117.6 1.2',
            'NOTE BUSTER_UTILANGLE 117.6 1.5 A|100:OP1.A A|100:P.A A|100:OP2.A'
        ]

        self.expected_lines_dist = [
            '# Restraint PO4==AA_0 dist dO1P4 1.487 0.010',
            'NOTE BUSTER_DISTANCE =1.487 0.020 A|100:OP1.A A|100:P.A'
        ]


class TuplePrinterSigmaOverrideTestCase(PrinterTestCase):

    def setUp(self):
        super(TuplePrinterSigmaOverrideTestCase, self).setUp()
        self.printer = TuplePrinter(override_sigma=True)

        self.expected_lines_angle = [
            ((1, 2, 3), 117.6, 3.0)
        ]

        self.expected_lines_dist = [
            ((1, 2), 1.487, 0.020)
        ]

    def test_get_angle(self):
        angle_text = self.printer.get_angle(self.restraints[0], self.restraints)

        for line in self.expected_lines_angle:
            self.assertEqual(line, angle_text)

    def test_get_dist(self):
        dist = self.printer.get_dist(self.restraints[1], self.restraints)

        for line in self.expected_lines_dist:
            self.assertEqual(line, dist)
