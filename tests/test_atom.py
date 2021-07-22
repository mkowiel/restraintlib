# -*- coding: utf-8 -*-
import math
from unittest import TestCase

from restraintlib.atom import Atom


class AtomTestCase(TestCase):
    def setUp(self):
        self.atom = Atom('A', '100', 'DT', 'C2', ' ', (1.0, 0.2, 0.3), 1)
        self.atom2 = Atom('A', '100', 'DT', 'C3', ' ', (1.0, 1.0, 1.0), 2)

    def test_str(self):
        expected = "chain: A res: 100 monomer: DT atom: C2 alt loc:  xyz: (1.0, 0.2, 0.3)"
        self.assertEqual(expected, str(self.atom))

    def test_cross(self):
        cross = Atom.cross(self.atom.atom_xyz, self.atom2.atom_xyz)
        self.assertEqual((0.2 - 0.3, 0.3 - 1.0, 1.0 - 0.2), cross)

    def test_sub(self):
        sub = Atom.sub(self.atom.atom_xyz, self.atom2.atom_xyz)
        self.assertEqual((0.0, -0.8, -0.7), sub)

    def test_dot(self):
        dot = Atom.dot(self.atom.atom_xyz, self.atom2.atom_xyz)
        self.assertEqual(1.5, dot)

    def test_lenght(self):
        length = Atom.lenght(self.atom.atom_xyz)
        self.assertAlmostEqual(math.sqrt(1.13), length, 5)

        length = Atom.lenght(self.atom2.atom_xyz)
        self.assertAlmostEqual(math.sqrt(3.0), length, 5)

    def test_mul_sca(self):
        vec = Atom.mul_sca(3, self.atom2.atom_xyz)
        self.assertEqual((3.0, 3.0, 3.0), vec)

    def test_normalize(self):
        vec = Atom.normalize(self.atom2.atom_xyz)
        self.assertEqual(1.0, Atom.lenght(vec))

    def test_det(self):
        det = Atom.det(self.atom.atom_xyz, self.atom2.atom_xyz, self.atom2.atom_xyz)
        self.assertAlmostEqual(0, det, 6)

        det = Atom.det(self.atom.atom_xyz, self.atom2.atom_xyz, (1, 1, 0.5))
        self.assertAlmostEqual(-0.4, det, 6)

    def test_dist(self):
        dist = self.atom.dist(self.atom2)
        self.assertAlmostEqual(math.sqrt(1.13), dist, 6)

    def test_angle(self):
        zero = Atom('', '1', '', '', '', (0.0, 0.0, 0.0), 1)
        one = Atom('', '2', '', '', '', (1.0, 0.0, 0.0), 2)
        diag = Atom('', '3', '', '', '', (1.0, 1.0, 0.0), 3)
        one_one_one = Atom('', '3', '', '', '', (1.0, 1.0, 1.0), 4)

        angle = one.angle(zero, diag)
        self.assertAlmostEqual(45, angle, 6)

        angle = diag.angle(zero, one_one_one)
        self.assertAlmostEqual(35.26, angle, 2)

    def test_torsion(self):
        a1 = Atom('', '1', '', '', '', (-1.0, -1.0, 0.0), 1)
        a2 = Atom('', '2', '', '', '', (-1.0, 0.0, 0.0), 2)
        a3 = Atom('', '3', '', '', '', (1.0, 0.0, 0.0), 3)
        a4 = Atom('', '3', '', '', '', (1.0, 1.0, 0.0), 4)

        torsion = a1.torsion(a2, a3, a4)
        self.assertAlmostEqual(180, torsion, 6)
