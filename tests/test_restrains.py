# -*- coding: utf-8 -*-
import math
from unittest import TestCase

from restraintlib.atom import Atom
from restraintlib.lib.PO4 import PO4_RESTRAINTS
from restraintlib.restraints import DistanceMeasure
from restraintlib.restraints import ConditionalRestraintItem
from restraintlib.restraints import ConditionItem
from restraintlib.restraints import ConditionalRestraint
from restraintlib.restraints import ConditionalRestraintList
from restraintlib.restraints import Restraint
from restraintlib.restraints import MonomerRestraintGroup


class DistanceMeasureTestCase(TestCase):

    def setUp(self):
        self.measure_euclidean = DistanceMeasure('euclidean', ['C1', 'C2'])
        self.measure_euclidean_angles = DistanceMeasure('euclidean_angles', ['C1', 'C2'])

    def test_euclidean(self):
        self.assertEqual(self.measure_euclidean.euclidean([0.0, 2.0], [0.0, 0.0]), 2.0)
        self.assertEqual(self.measure_euclidean.euclidean([-2.0, 2.0], [0.0, 0.0]), 2.0*math.sqrt(2.0))

    def test_euclidean_angles(self):
        self.assertEqual(self.measure_euclidean_angles.euclidean_angles([0.0], [360.0]), 0.0)
        self.assertEqual(self.measure_euclidean_angles.euclidean_angles([0.0, 2.0], [360.0, 0.0]), 2.0)
        self.assertEqual(
            self.measure_euclidean_angles.euclidean_angles([-2.0, 2.0], [360.0, 0.0]), 2.0 * math.sqrt(2.0)
        )

    def test_distance(self):
        pass


class ConditionalRestraintItemTestCase(TestCase):

    def setUp(self):
        self.item = ConditionalRestraintItem(None, None, None, None, None, value_dist=None, sigma_dist=None)


class ConditionItemTestCase(TestCase):
    def setUp(self):
        self.item_pseudorotation = ConditionItem(
            'pseudorotation',
            'test_pseudorotation',
            ["C1'", "C2'", "C3'", "C4'", "O4'"],
            162,
            18,
        )

    def test_fix_torsion(self):
        self.assertEqual(self.item_pseudorotation.fix_torsion(0), 0)
        self.assertEqual(self.item_pseudorotation.fix_torsion(40), 40)
        self.assertEqual(self.item_pseudorotation.fix_torsion(190), -170)
        self.assertEqual(self.item_pseudorotation.fix_torsion(-190), 170)
        self.assertEqual(self.item_pseudorotation.fix_torsion(180), 180)
        self.assertEqual(self.item_pseudorotation.fix_torsion(-180), -180)
        self.assertEqual(self.item_pseudorotation.fix_torsion(-350.5), 9.5)
        self.assertEqual(self.item_pseudorotation.fix_torsion(-370.5), -10.5)
        self.assertEqual(self.item_pseudorotation.fix_torsion(-370.5-180), -10.5+180)
        self.assertEqual(self.item_pseudorotation.fix_torsion(350.5), -9.5)
        self.assertEqual(self.item_pseudorotation.fix_torsion(370.5), 10.5)
        self.assertEqual(self.item_pseudorotation.fix_torsion(180+370.5), 10.5-180)

    def test__pseudorotation(self):
        thetas = [
            [-7.33, -8.64, 19.28, -24.05, 19.83],
            [-34.12, 40.64, -32.11, 12.42, 12.59],
            [-25.11, 35.39, -31.71, 18.12, 4.33],
        ]
        expected_ps = [35.66, 142.83, 154.42]
        expected_tms = [23.73, 40.30, 35.16]

        for theta, expected_p, expected_tm in zip(thetas, expected_ps, expected_tms):
            p, tm = self.item_pseudorotation._pseudorotation(theta)
            self.assertAlmostEqual(p, expected_p, delta=0.01)
            self.assertAlmostEqual(tm, expected_tm, delta=0.01)

    def test__pseudorotation_with_sd(self):
        thetas = [
            [-7.33, -8.64, 19.28, -24.05, 19.83],
            [-34.12, 40.64, -32.11, 12.42, 12.59],
            [-25.11, 35.39, -31.71, 18.12, 4.33],
        ]
        expected_ps = [35.66, 142.83, 154.42]
        expected_tms = [23.73, 40.30, 35.16]

        for theta, expected_p, expected_tm in zip(thetas, expected_ps, expected_tms):
            p, sd_p, tm, sd_tm = self.item_pseudorotation._pseudorotation_with_sd(theta)
            self.assertAlmostEqual(p, expected_p, delta=0.6)
            self.assertAlmostEqual(tm, expected_tm, delta=0.9)


class ConditionalRestraintTestCase(TestCase):

    def setUp(self):
        self.restraint = ConditionalRestraint(None, None, None)


class ConditionalRestraintListTestCase(TestCase):

    def setUp(self):
        self.restraint_list = ConditionalRestraintList(PO4_RESTRAINTS[4:6])

    def test_find_restraint_closest(self):
        atoms = {
            "OP1": Atom('A', 20, 'C', "OP1", '', [1, 0, 0], 1),
            "OP2": Atom('A', 20, 'C', "OP2", '', [0, 1, 0], 2),
            "P":   Atom('A', 20, 'C', "P",   '', [0, 0, 0], 3),
            "O3'": Atom('A', 20, 'C', "O3'", '', [1, 0, 1], 4),
            "O5'": Atom('A', 20, 'C', "O5'", '', [0, 1, 1.1], 5),
        }
        distance_measure = DistanceMeasure('euclidean_angles', ['aO1O2', 'aO1O3', 'aO1O5', 'aO2O3', 'aO2O5', 'aO3O5'])
        closest = self.restraint_list.find_restraint_closest(atoms, distance_measure)
        self.assertEqual(closest.name, 'PO4==AS_0')


class RestraintTestCase(TestCase):

    def setUp(self):
        self.restraint = Restraint(None, None, None, None, None, None, None, None)


class MonomerRestraintGroupTestCase(TestCase):

    def setUp(self):
        distance_measure = {
            'measure': 'euclidean_angles',
            'restraint_names': []
        }

        restraints = [
            {
                "conditions": [
                    ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], -66.636, 7.779],
                    ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], 171.37, 14.971]
                ],
                "name": "PO4==AA_0",
                "restraints": [
                    ["angle", "aO1O2", ["OP1", "P", "OP2"], 117.6, 1.2],
                    ["angle", "aO1O3", ["OP1", "P", "O3'"], 106.2, 1.1],
                    ["angle", "aO1O5", ["OP1", "P", "O5'"], 110.2, 1.3],
                    ["angle", "aO2O3", ["OP2", "P", "O3'"], 112.2, 1.0],
                    ["angle", "aO2O5", ["OP2", "P", "O5'"], 109.3, 0.9],
                    ["angle", "aO3O5", ["O3'", "P", "O5'"], 99.9, 0.7],
                    ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.2, 1.5],
                    ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 121.7, 3.0],
                    ["dist", "dO1P4", ["OP1", "P"], 1.487, 0.01],
                    ["dist", "dO2P4", ["OP2", "P"], 1.483, 0.01],
                    ["dist", "dO3P4", ["O3'", "P"], 1.601, 0.008],
                    ["dist", "dO5P4", ["O5'", "P"], 1.591, 0.004],
                    ["dist", "dO3C3", ["O3'", "C3'"], 1.422, 0.010],
                    ["dist", "dO5C5", ["O5'", "C5'"], 1.428, 0.013]
                ]
            },
        ]

        self.restraint = MonomerRestraintGroup(None, None, {}, None, None, None, restraints, distance_measure, distance_measure)

    def test_register_res_id(self):
        self.restraint.register_res_id('A', '1')
        self.restraint.register_res_id('A', '2')
        self.restraint.register_res_id('A', '3')

        self.assertEqual(self.restraint.prev_res_id('A', '1'), None)
        self.assertEqual(self.restraint.next_res_id('A', '1'), '2')

        self.assertEqual(self.restraint.prev_res_id('A', '2'), '1')
        self.assertEqual(self.restraint.next_res_id('A', '2'), '3')

        self.assertEqual(self.restraint.prev_res_id('A', '3'), '2')
        self.assertEqual(self.restraint.next_res_id('A', '3'), None)
