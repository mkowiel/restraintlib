# coding: utf-8
from __future__ import print_function

import math

__author__ = "Marcin Kowiel, Dariusz Brzezinski"


class Atom(object):
    # use slots to minimize memory footprint
    __slots__ = (
        'chain_id', 'res_id', 'res_name', 'atom_name', 'alt_loc', 'atom_xyz', 'i_seq'
    )

    def __init__(self, chain_id, res_id, res_name, atom_name, alt_loc, atom_xyz, i_seq):
        """
        Represents an atom from a pdb file, for which restraints (e.g., bond length or angle) will be created.
        :param chain_id: chain id
        :param res_id: residual id
        :param res_name: residual name
        :param atom_name: atom name
        :param alt_loc: alternate location
        :param atom_xyz: xyz coordinates of the atom
        """
        self.chain_id = chain_id.strip()
        self.res_id = res_id
        self.res_name = res_name.strip()
        self.atom_name = atom_name.strip()
        self.alt_loc = alt_loc.strip()
        self.atom_xyz = atom_xyz
        self.i_seq = i_seq

    def __str__(self):
        """
        String representation of the atom.
        """
        return "chain: {} res: {} monomer: {} atom: {} alt loc: {} xyz: {}".format(
            self.chain_id, self.res_id, self.res_name, self.atom_name, self.alt_loc, self.atom_xyz)

    def __eq__(self, other):
        return (
            self.chain_id == other.chain_id and
            self.res_id == other.res_id and
            self.res_name == other.res_name and
            self.atom_name == other.atom_name and
            self.alt_loc == other.alt_loc and
            self.atom_xyz == other.atom_xyz
        )

    @classmethod
    def cross3(cls, x1, y1, z1, x2, y2, z2):
        return (y1*z2-z1*y2, -x1*z2+z1*x2, x1*y2-y1*x2)

    @classmethod
    def cross(cls, a, b):
        return cls.cross3(a[0], a[1], a[2], b[0], b[1], b[2])

    @classmethod
    def sub(cls, a, b):
        return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

    @classmethod
    def dot(cls, a, b):
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

    @classmethod
    def lenght(cls, a):
        return math.sqrt(cls.dot(a, a))

    @classmethod
    def mul_sca(cls, sca, vec):
        return (sca*vec[0],sca*vec[1], sca*vec[2])

    @classmethod
    def normalize(cls, vec):
        vec_lenght = cls.lenght(vec)
        return cls.mul_sca(1.0/vec_lenght, vec)

    @classmethod
    def det(cls, a, b, c):
        return cls.dot(cls.cross(a, b), c)

    def dist(self, atom):
        dx = self.atom_xyz[0]-atom.atom_xyz[0]
        dy = self.atom_xyz[1]-atom.atom_xyz[1]
        dz = self.atom_xyz[2]-atom.atom_xyz[2]
        return math.sqrt(dx*dx+dy*dy+dz*dz)

    def angle(self, atom2, atom3):
        a = self.sub(self.atom_xyz, atom2.atom_xyz)
        b = self.sub(atom3.atom_xyz, atom2.atom_xyz)
        c = self.dot(a, b)
        d = self.lenght(a)*self.lenght(b)
        angle = math.acos(c/d)
        return math.degrees(angle)

    def torsion(self, atom2, atom3, atom4):
        """
        Return torsion of self -> atom2, atom2 -> atom3, atom3 -> atom4 vectors.
        :param atom2:
        :param atom3:
        :param atom4:
        :return: torsion angle in degrees
        """
        vector1 = self.sub(atom2.atom_xyz, self.atom_xyz)
        vector2 = self.sub(atom3.atom_xyz, atom2.atom_xyz)
        vector3 = self.sub(atom4.atom_xyz, atom3.atom_xyz)

        #torsion angle = atan2 (|b2|b1 . (b2 x b3), (b1 x b2) . (b2 x b3)
        #source: http://en.wikipedia.org/wiki/Dihedral_angle

        y2 = self.cross(vector2, vector3)
        y = self.lenght(vector2)*self.dot(vector1, y2)

        x1 = self.cross(vector1, vector2)
        x = self.dot(x1, y2)

        torsion = math.atan2(y, x)
        torsion = math.degrees(torsion)
        return torsion