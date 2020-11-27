# -*- coding: utf-8 -*-
import math
from unittest import TestCase

from restraintlib.atom import Atom
from restraintlib.lib.PO4 import PO4_RESTRAINTS
from restraintlib.restrains import DistanceMeasure
from restraintlib.restrains import ConditionalRestraintItem
from restraintlib.restrains import ConditionItem
from restraintlib.restrains import ConditionalRestraint
from restraintlib.restrains import ConditionalRestraintList
from restraintlib.restrains import Restraint
from restraintlib.restrains import MonomerRestraintGroup


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
            "OP1": Atom('A', 20, 'C', "OP1", '', [1, 0, 0]),
            "OP2": Atom('A', 20, 'C', "OP2", '', [0, 1, 0]),
            "P":   Atom('A', 20, 'C', "P",   '', [0, 0, 0]),
            "O3'": Atom('A', 20, 'C', "O3'", '', [1, 0, 1]),
            "O5'": Atom('A', 20, 'C', "O5'", '', [0, 1, 1.1]),
        }
        distance_measure = DistanceMeasure('euclidean_angles', ['aO1O2', 'aO1O3', 'aO1O5', 'aO2O3', 'aO2O5', 'aO3O5'])
        closest = self.restraint_list.find_restraint_closest(atoms, distance_measure)
        self.assertEqual(closest.name, 'PO4==AS_0')


class RestraintTestCase(TestCase):

    def setUp(self):
        self.restraint = Restraint(None, None, None, None, None, None, None, None)


class MonomerRestraintGroupTestCase(TestCase):

    def setUp(self):
        self.restraint = MonomerRestraintGroup(None, None, None, None, None, None, None, None)
