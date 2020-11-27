# -*- coding: utf-8 -*-
import os
import re
from StringIO import StringIO

from unittest import TestCase
from unittest import expectedFailure

from restraintlib.restrainlib import AllowedRestraintsConfig
from restraintlib.restrainlib import RestrainLib
from restraintlib.restrains import VERSION


class RestrainLibTestCase(TestCase):
    def setUp(self):
        self.buffer = StringIO()
        self.lib = RestrainLib(log_stream=self.buffer)
        self.pdb_3p4j = os.path.join(os.path.dirname(os.path.abspath(__file__)), '3p4j.pdb')
        self.pdb_1d8g = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1d8g.pdb')
        self.pdb_ig = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'IG.pdb')
        self.pdb_4r15 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '4r15.pdb')
        self.pdb_disorder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'disorder.pdb')
        self.mmcif_3p4j = os.path.join(os.path.dirname(os.path.abspath(__file__)), '3p4j.cif')
        self.mmcif_1d8g = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1d8g.cif')
        self.mmcif_3ssf = os.path.join(os.path.dirname(os.path.abspath(__file__)), '3ssf.cif')
        self.restraints_config = AllowedRestraintsConfig()

    def assert_citing(self, data):
        self.assertIn(VERSION, data)
        self.assertIn("Nucleic Acids Res.", data)
        self.assertIn("10.1093/nar/gkw717", data)

    def assert_adenine(self, data):
        # adenine
        self.assertIn("dN1C2", data)
        self.assertIn("dC2N3", data)
        self.assertIn("dN3C4", data)
        self.assertIn("dC4C5", data)
        self.assertIn("dC5C6", data)
        self.assertIn("dC6N1", data)
        self.assertIn("dC5N7", data)
        self.assertIn("dN7C8", data)
        self.assertIn("dC8N9", data)
        self.assertIn("dN9C4", data)
        self.assertIn("dC6N6", data)

        self.assertIn("aC6N1C2", data)
        self.assertIn("aN1C2N3", data)
        self.assertIn("aC2N3C4", data)
        self.assertIn("aN3C4C5", data)
        self.assertIn("aC4C5C6", data)
        self.assertIn("aC5C6N1", data)
        self.assertIn("aN3C4N9", data)
        self.assertIn("aC6C5N7", data)
        self.assertIn("aC5C4N9", data)
        self.assertIn("aC4N9C8", data)
        self.assertIn("aN9C8N7", data)
        self.assertIn("aC8N7C5", data)
        self.assertIn("aN7C5C4", data)
        self.assertIn("aN6C6N1", data)
        self.assertIn("aN6C6C5", data)

    def assert_guanine(self, data):
        self.assertIn("dN1C2", data)
        self.assertIn("dC2N3", data)
        self.assertIn("dN3C4", data)
        self.assertIn("dC4C5", data)
        self.assertIn("dC5C6", data)
        self.assertIn("dC6N1", data)
        self.assertIn("dC5N7", data)
        self.assertIn("dN7C8", data)
        self.assertIn("dC8N9", data)
        self.assertIn("dN9C4", data)
        self.assertIn("dC6O6", data)
        self.assertIn("dC2N2", data)

        self.assertIn("aC6N1C2", data)
        self.assertIn("aN1C2N3", data)
        self.assertIn("aC2N3C4", data)
        self.assertIn("aN3C4C5", data)
        self.assertIn("aC4C5C6", data)
        self.assertIn("aC5C6N1", data)
        self.assertIn("aN3C4N9", data)
        self.assertIn("aC6C5N7", data)
        self.assertIn("aC5C4N9", data)
        self.assertIn("aC4N9C8", data)
        self.assertIn("aN9C8N7", data)
        self.assertIn("aC8N7C5", data)
        self.assertIn("aN7C5C4", data)
        self.assertIn("aO6C6N1", data)
        self.assertIn("aO6C6C5", data)
        self.assertIn("aN2C2N1", data)
        self.assertIn("aN2C2N3", data)

    def assert_uracil(self, data):
        self.assertIn("dN1C2", data)
        self.assertIn("dC2N3", data)
        self.assertIn("dN3C4", data)
        self.assertIn("dC4C5", data)
        self.assertIn("dC5C6", data)
        self.assertIn("dC6N1", data)
        self.assertIn("dC2O2", data)
        self.assertIn("dC4O4", data)

        self.assertIn("aC6N1C2", data)
        self.assertIn("aN1C2N3", data)
        self.assertIn("aC2N3C4", data)
        self.assertIn("aN3C4C5", data)
        self.assertIn("aC4C5C6", data)
        self.assertIn("aC5C6N1", data)
        self.assertIn("aO2C2N1", data)
        self.assertIn("aO2C2N3", data)
        self.assertIn("aO4C4C5", data)
        self.assertIn("aO4C4N3", data)

        self.assertIn("1.381", data)
        self.assertIn("1.373", data)
        self.assertIn("1.381", data)
        self.assertIn("1.432", data)
        self.assertIn("1.337", data)
        self.assertIn("1.374", data)
        self.assertIn("1.219", data)
        self.assertIn("1.231", data)

        self.assertIn("121.1", data)
        self.assertIn("114.9", data)
        self.assertIn("127.0", data)
        self.assertIn("114.5", data)
        self.assertIn("119.7", data)
        self.assertIn("122.7", data)
        self.assertIn("122.8", data)
        self.assertIn("122.3", data)
        self.assertIn("126.0", data)
        self.assertIn("119.5", data)

    def assert_thymine(self, data):
        self.assertIn("dN1C2", data)
        self.assertIn("dC2N3", data)
        self.assertIn("dN3C4", data)
        self.assertIn("dC4C5", data)
        self.assertIn("dC5C6", data)
        self.assertIn("dC6N1", data)
        self.assertIn("dC2O2", data)
        self.assertIn("dC4O4", data)
        self.assertIn("dC7C5", data)

        self.assertIn("aC6N1C2", data)
        self.assertIn("aN1C2N3", data)
        self.assertIn("aC2N3C4", data)
        self.assertIn("aN3C4C5", data)
        self.assertIn("aC4C5C6", data)
        self.assertIn("aC5C6N1", data)
        self.assertIn("aO2C2N1", data)
        self.assertIn("aO2C2N3", data)
        self.assertIn("aO4C4C5", data)
        self.assertIn("aO4C4N3", data)
        self.assertIn("aC7C5C4", data)
        self.assertIn("aC7C5C6", data)

        self.assertIn("1.376", data)
        self.assertIn("1.372", data)
        self.assertIn("1.382", data)
        self.assertIn("1.446", data)
        self.assertIn("1.340", data)
        self.assertIn("1.381", data)
        self.assertIn("1.222", data)
        self.assertIn("1.229", data)
        self.assertIn("1.498", data)

        self.assertIn("121.2", data)
        self.assertIn("114.7", data)
        self.assertIn("127.1", data)
        self.assertIn("115.2", data)
        self.assertIn("118.1", data)
        self.assertIn("123.6", data)
        self.assertIn("123.0", data)
        self.assertIn("122.3", data)
        self.assertIn("125.0", data)
        self.assertIn("119.8", data)
        self.assertIn("118.7", data)
        self.assertIn("123.2", data)

    def assert_cytosine(self, data):
        self.assertIn("dN1C2", data)
        self.assertIn("dC2N3", data)
        self.assertIn("dN3C4", data)
        self.assertIn("dC4C5", data)
        self.assertIn("dC5C6", data)
        self.assertIn("dC6N1", data)
        self.assertIn("dC2O2", data)
        self.assertIn("dC4N4", data)

        self.assertIn("aC6N1C2", data)
        self.assertIn("aN1C2N3", data)
        self.assertIn("aC2N3C4", data)
        self.assertIn("aN3C4C5", data)
        self.assertIn("aC4C5C6", data)
        self.assertIn("aC5C6N1", data)
        self.assertIn("aO2C2N1", data)
        self.assertIn("aO2C2N3", data)
        self.assertIn("aN4C4C5", data)
        self.assertIn("aN4C4N3", data)

        self.assertIn("1.395", data)
        self.assertIn("1.353", data)
        self.assertIn("1.337", data)
        self.assertIn("1.424", data)
        self.assertIn("1.338", data)
        self.assertIn("1.365", data)
        self.assertIn("1.240", data)
        self.assertIn("1.330", data)

        self.assertIn("120.3", data)
        self.assertIn("119.1", data)
        self.assertIn("120.1", data)
        self.assertIn("121.6", data)
        self.assertIn("117.5", data)
        self.assertIn("121.2", data)
        self.assertIn("118.8", data)
        self.assertIn("122.0", data)
        self.assertIn("120.3", data)
        self.assertIn("118.1", data)

    def assert_isocytosine(self, data):
        self.assertIn("dN1C2", data)
        self.assertIn("dC2N3", data)
        self.assertIn("dN3C4", data)
        self.assertIn("dC4C5", data)
        self.assertIn("dC5C6", data)
        self.assertIn("dC6N1", data)
        self.assertIn("dC2N2", data)
        self.assertIn("dC4O4", data)

        self.assertIn("aC6N1C2", data)
        self.assertIn("aN1C2N3", data)
        self.assertIn("aC2N3C4", data)
        self.assertIn("aN3C4C5", data)
        self.assertIn("aC4C5C6", data)
        self.assertIn("aC5C6N1", data)
        self.assertIn("aN2C2N1", data)
        self.assertIn("aN2C2N3", data)
        self.assertIn("aO4C4C5", data)
        self.assertIn("aO4C4N3", data)

        self.assertIn("1.385", data)
        self.assertIn("1.324", data)
        self.assertIn("1.375", data)
        self.assertIn("1.461", data)
        self.assertIn("1.350", data)
        self.assertIn("1.382", data)
        self.assertIn("1.340", data)
        self.assertIn("1.244", data)

        self.assertIn("118.2", data)
        self.assertIn("122.1", data)
        self.assertIn("121.7", data)
        self.assertIn("117.3", data)
        self.assertIn("119.0", data)
        self.assertIn("121.7", data)
        self.assertIn("118.8", data)
        self.assertIn("119.0", data)
        self.assertIn("121.9", data)
        self.assertIn("120.8", data)

    def assert_isoguanine(self, data):
        self.assertIn("dN1C2", data)
        self.assertIn("dC2N3", data)
        self.assertIn("dN3C4", data)
        self.assertIn("dC4C5", data)
        self.assertIn("dC5C6", data)
        self.assertIn("dC6N1", data)
        self.assertIn("dC5N7", data)
        self.assertIn("dN7C8", data)
        self.assertIn("dC8N9", data)
        self.assertIn("dN9C4", data)
        self.assertIn("dC6N6", data)
        self.assertIn("dC2O2", data)

        self.assertIn("aC6N1C2", data)
        self.assertIn("aN1C2N3", data)
        self.assertIn("aC2N3C4", data)
        self.assertIn("aN3C4C5", data)
        self.assertIn("aC4C5C6", data)
        self.assertIn("aC5C6N1", data)
        self.assertIn("aN3C4N9", data)
        self.assertIn("aC6C5N7", data)
        self.assertIn("aC5C4N9", data)
        self.assertIn("aC4N9C8", data)
        self.assertIn("aN9C8N7", data)
        self.assertIn("aC8N7C5", data)
        self.assertIn("aN7C5C4", data)
        self.assertIn("aN6C6N1", data)
        self.assertIn("aN6C6C5", data)
        self.assertIn("aO2C2N1", data)
        self.assertIn("aO2C2N3", data)

        self.assertIn("1.419", data)
        self.assertIn("1.351", data)
        self.assertIn("1.334", data)
        self.assertIn("1.407", data)
        self.assertIn("1.401", data)
        self.assertIn("1.361", data)
        self.assertIn("1.388", data)
        self.assertIn("1.312", data)
        self.assertIn("1.390", data)
        self.assertIn("1.378", data)
        self.assertIn("1.337", data)
        self.assertIn("1.254", data)

        self.assertIn("126.0", data)
        self.assertIn("119.6", data)
        self.assertIn("114.0", data)
        self.assertIn("129.3", data)
        self.assertIn("116.5", data)
        self.assertIn("114.5", data)
        self.assertIn("126.3", data)
        self.assertIn("132.0", data)
        self.assertIn("104.4", data)
        self.assertIn("106.7", data)
        self.assertIn("113.4", data)
        self.assertIn("104.1", data)
        self.assertIn("111.4", data)
        self.assertIn("119.8", data)
        self.assertIn("125.7", data)
        self.assertIn("116.6", data)
        self.assertIn("123.8", data)

    def _assert_sugar(self, data):
        self.assertIn("dC1'C2'", data)
        self.assertIn("dC2'C3'", data)
        self.assertIn("dC3'C4'", data)
        self.assertIn("dC4'O4'", data)
        self.assertIn("dC1'O4'", data)
        self.assertIn("dC4'C5'", data)
        self.assertIn("dC1'N1", data)
        self.assertIn("dC3'O3'", data)
        self.assertIn("dC5'O5'", data)

        self.assertIn("aN1C1'C2'", data)
        self.assertIn("aN9C1'C2'", data)
        self.assertIn("aN1C1'O4'", data)
        self.assertIn("aN9C1'O4'", data)
        self.assertIn("aC1'N1C2", data)
        self.assertIn("aC1'N9C4", data)
        self.assertIn("aC1'N1C6", data)
        self.assertIn("aC1'N9C8", data)
        self.assertIn("aC1'C2'C3'", data)
        self.assertIn("aC2'C3'C4'", data)
        self.assertIn("aC3'C4'O4'", data)
        self.assertIn("aC1'O4'C4'", data)
        self.assertIn("aC4'C3'O3'", data)
        self.assertIn("aC2'C3'O3'", data)
        self.assertIn("aC2'C1'O4'", data)
        self.assertIn("aC3'C4'C5'", data)
        self.assertIn("aC5'C4'O4'", data)
        self.assertIn("aC4'C5'O5'", data)

    def assert_deoxyribose(self, data):
        self._assert_sugar(data)
        self.assertNotIn("dC2'O2'", data)
        self.assertNotIn("aC1'C2'O2'", data)
        self.assertNotIn("aC3'C2'O2'", data)

    def assert_ribose(self, data):
        self._assert_sugar(data)
        self.assertIn("dC2'O2'", data)
        self.assertIn("aC1'C2'O2'", data)
        self.assertIn("aC3'C2'O2'", data)

    def assert_3p4j(self, data):
        self.assert_citing(data)

        self.assertIn("AA_3", data)
        self.assertIn("AS_3", data)
        self.assert_guanine(data)
        self.assert_cytosine(data)
        self.assert_deoxyribose(data)

    def assert_1d8g(self, data):
        self.assert_citing(data)

        self.assertIn("AA_1", data)
        self.assertIn("AS_0", data)
        self.assertIn("AS_1", data)
        self.assert_adenine(data)
        self.assert_guanine(data)
        self.assert_thymine(data)
        self.assert_cytosine(data)
        self.assert_deoxyribose(data)

    def assert_3ssf(self, data):
        self.assert_citing(data)

        self.assert_adenine(data)
        self.assert_guanine(data)
        self.assert_thymine(data)
        self.assert_cytosine(data)
        self.assert_ribose(data)

    def assert_iso2(self, data):
        self.assert_citing(data)

        self.assert_adenine(data)
        self.assert_guanine(data)
        self.assert_cytosine(data)
        self.assert_isocytosine(data)
        self.assert_isoguanine(data)
        self.assert_ribose(data)

    def assert_4r15(self, data):
        self.assert_citing(data)

        self.assertNotIn("dN1C2", data)
        self.assertNotIn("dC2N3", data)
        self.assertNotIn("dN3C4", data)
        self.assertNotIn("dC4C5", data)
        self.assertNotIn("dC5C6", data)
        self.assertNotIn("dC6N1", data)
        self.assertNotIn("dC5N7", data)
        self.assertNotIn("dN7C8", data)
        self.assertNotIn("dC8N9", data)
        self.assertNotIn("dN9C4", data)
        self.assertNotIn("dC6O6", data)
        self.assertNotIn("dC2N2", data)

        self.assertNotIn("aC6N1C2", data)
        self.assertNotIn("aN1C2N3", data)
        self.assertNotIn("aC2N3C4", data)
        self.assertNotIn("aN3C4C5", data)
        self.assertNotIn("aC4C5C6", data)
        self.assertNotIn("aC5C6N1", data)
        self.assertNotIn("aN3C4N9", data)
        self.assertNotIn("aC6C5N7", data)
        self.assertNotIn("aC5C4N9", data)
        self.assertNotIn("aC4N9C8", data)
        self.assertNotIn("aN9C8N7", data)
        self.assertNotIn("aC8N7C5", data)
        self.assertNotIn("aN7C5C4", data)
        self.assertNotIn("aO6C6N1", data)
        self.assertNotIn("aO6C6C5", data)
        self.assertNotIn("aN2C2N1", data)
        self.assertNotIn("aN2C2N3", data)

        self.assertNotIn("dN1C2", data)
        self.assertNotIn("dC2N3", data)
        self.assertNotIn("dN3C4", data)
        self.assertNotIn("dC4C5", data)
        self.assertNotIn("dC5C6", data)
        self.assertNotIn("dC6N1", data)
        self.assertNotIn("dC2O2", data)
        self.assertNotIn("dC4N4", data)

        self.assertNotIn("aC6N1C2", data)
        self.assertNotIn("aN1C2N3", data)
        self.assertNotIn("aC2N3C4", data)
        self.assertNotIn("aN3C4C5", data)
        self.assertNotIn("aC4C5C6", data)
        self.assertNotIn("aC5C6N1", data)
        self.assertNotIn("aO2C2N1", data)
        self.assertNotIn("aO2C2N3", data)
        self.assertNotIn("aN4C4C5", data)
        self.assertNotIn("aN4C4N3", data)

        self.assert_deoxyribose(data)

    def assert_length_for_shelx(self, data, max_line_size=80):
        for line in data.split('\n'):
            self.assertLessEqual(len(line), max_line_size, msg="Line too long {}".format(line))

    def test_produce_restrains_phenix_3p4j_pdb(self):
        self.lib.produce_restrains(self.pdb_3p4j, 'P', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_3p4j(data)

    def test_produce_restrains_shelx_3p4j_pdb(self):
        self.lib.produce_restrains(self.pdb_3p4j, 'S', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_3p4j(data)
        self.assert_length_for_shelx(data)

    def test_produce_restrains_refmac_3p4j_pdb(self):
        self.lib.produce_restrains(self.pdb_3p4j, 'R', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_3p4j(data)

    def test_produce_restrains_phenix_3p4j_mmcif(self):
        self.lib.produce_restrains(self.mmcif_3p4j, 'P', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_3p4j(data)

    def test_produce_restrains_shelx_3p4j_mmcif(self):
        self.lib.produce_restrains(self.mmcif_3p4j, 'S', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_3p4j(data)
        self.assert_length_for_shelx(data)

    def test_produce_restrains_refmac_3p4j_mmcif(self):
        self.lib.produce_restrains(self.mmcif_3p4j, 'R', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_3p4j(data)

    def test_produce_restrains_phenix_1d8g_pdb(self):
        self.lib.produce_restrains(self.pdb_1d8g, 'P', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_1d8g(data)

    def test_produce_restrains_shelx_1d8g_pdb(self):
        self.lib.produce_restrains(self.pdb_1d8g, 'S', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_1d8g(data)
        self.assert_length_for_shelx(data)

    def test_produce_restrains_refmac_1d8g_pdb(self):
        self.lib.produce_restrains(self.pdb_1d8g, 'R', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_1d8g(data)

    def test_produce_restrains_phenix_3ssf_mmcif(self):
        self.lib.produce_restrains(self.mmcif_3ssf, 'P', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_3ssf(data)

    def test_produce_restrains_phenix_3p4j_pdb_mmcif_same(self):
        self.lib.produce_restrains(self.pdb_3p4j, 'P', self.restraints_config)
        data_pdb = self.buffer.getvalue()
        self.buffer = StringIO()
        self.lib = RestrainLib(log_stream=self.buffer)
        self.lib.produce_restrains(self.mmcif_3p4j, 'P', self.restraints_config)
        data_cif = self.buffer.getvalue()
        data_pdb = [line for line in data_pdb.splitlines() if not line.startswith("#")]
        data_cif = [line for line in data_cif.splitlines() if not line.startswith("#")]
        self.assertEqual(data_pdb, data_cif)

    def test_produce_restrains_refmac_iso2_pdb(self):
        self.lib.produce_restrains(self.pdb_ig, 'R', self.restraints_config)
        data = self.buffer.getvalue()
        self.assert_iso2(data)
        self.assertIn("exte dist first chain D resi 11 atom C1' second chain D resi 11 atom N9 value 1.", data)
        self.assertIn("ribose_purine_terminal_C5", data)
        self.assertIn("ribose_pyrimidine_terminal_C3", data)

    # fail since there are multiple outputs for surrounding of disorder atoms
    @expectedFailure
    def test_produce_restrains_phenix_disorder_pdb(self):
        self.lib.produce_restrains(self.pdb_disorder, 'P', self.restraints_config)
        data = self.buffer.getvalue()
        # if there are no duplicates there should be 2
        self.assertEqual(
            4, len([m.start() for m in re.finditer('atom_selection_1 = chain B and resid 20 and name OP1', data)])
        )

    def test_produce_restrains_shelx_4r15_pdb(self):
        restraints_config = AllowedRestraintsConfig(
            po4=False,
            po4terminal=False,
            nucleicacidbases=False,
            nucleicacidisobases=False,
        )
        self.lib.produce_restrains(self.pdb_4r15, 'S', restraints_config)
        data = self.buffer.getvalue()
        self.assert_4r15(data)
        self.assert_length_for_shelx(data)
