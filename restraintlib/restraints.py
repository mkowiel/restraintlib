# to be run with cctbx.python
import iotbx.pdb
import math
import os
import sys
from collections import defaultdict
from sklearn.externals import joblib

from atom import Atom
from lib import PO4
from lib import PO4_terminal_C3
from lib import PO4_terminal_C5
from lib import deoxyribose_purine
from lib import deoxyribose_purine_terminal_C3
from lib import deoxyribose_purine_terminal_C5
from lib import deoxyribose_pyrimidine
from lib import deoxyribose_pyrimidine_terminal_C3
from lib import deoxyribose_pyrimidine_terminal_C5
from lib import nucleic_acid_bases
from lib import nucleic_acid_isobases
from lib import ribose_purine
from lib import ribose_purine_terminal_C3
from lib import ribose_purine_terminal_C5
from lib import ribose_pyrimidine
from lib import ribose_pyrimidine_terminal_C3
from lib import ribose_pyrimidine_terminal_C5
from printer import CsvPrinter
from printer import PhenixPrinter
from printer import RefmacPrinter
from printer import ShelxPrinter

from cctbx.array_family import flex

VERSION = '2019.10.5'


class DistanceMeasure(object):
    # use slots to minimize memory footprint
    __slots__ = (
        'measure', 'restraint_names',
    )

    def __init__(self, measure, restraint_names):
        self.measure = getattr(self, measure)
        self.restraint_names = restraint_names

    @classmethod
    def euclidean(cls, vector1, vector2):
        if len(vector1) == 0 or len(vector2) == 0:
            return None
        if len(vector1) != len(vector2):
            raise Exception('uneven number of elements')

        dist_sq_sum = 0.0
        for a, b in zip(vector1, vector2):
            diff = a-b
            dist_sq_sum += diff*diff
        return math.sqrt(dist_sq_sum)

    @classmethod
    def euclidean_angles(cls, vector1, vector2):
        if len(vector1) == 0 or len(vector2) == 0:
            return None
        if len(vector1) != len(vector2):
            raise Exception('uneven number of elements')

        dist_sq_sum = 0.0
        for a, b in zip(vector1, vector2):
            diff = a-b
            diff = ConditionItem.fix_torsion(diff)
            dist_sq_sum += diff*diff
        return math.sqrt(dist_sq_sum)

    @classmethod
    def atoms_dist(cls, restraint, atoms):
        atom0 = atoms[restraint.atom_names[0]]
        atom1 = atoms[restraint.atom_names[1]]
        return atom0.dist(atom1)

    @classmethod
    def atoms_angle(cls, restraint, atoms):
        atom0 = atoms[restraint.atom_names[0]]
        atom1 = atoms[restraint.atom_names[1]]
        atom2 = atoms[restraint.atom_names[2]]
        return atom0.angle(atom1, atom2)

    @classmethod
    def atoms_torsion(cls, restraint, atoms):
        atom0 = atoms[restraint.atom_names[0]]
        atom1 = atoms[restraint.atom_names[1]]
        atom2 = atoms[restraint.atom_names[2]]
        atom3 = atoms[restraint.atom_names[3]]

        torsion = atom0.torsion(atom1, atom2, atom3)
        torsion = ConditionItem.fix_torsion(torsion)
        # TODO
        # even though it is normalized it is quite far from -179 to 179, but the distance should be taken mod 360
        return torsion

    @classmethod
    def atoms_pseudorotation(cls, restraint, atoms):
        C1 = atoms[restraint.atom_names[0]]
        C2 = atoms[restraint.atom_names[1]]
        C3 = atoms[restraint.atom_names[2]]
        C4 = atoms[restraint.atom_names[3]]
        O4 = atoms[restraint.atom_names[4]]

        theta = [0.0, 0.0, 0.0, 0.0, 0.0]
        theta[0] = ConditionItem.fix_torsion(C4.torsion(O4, C1, C2))
        theta[1] = ConditionItem.fix_torsion(O4.torsion(C1, C2, C3))
        theta[2] = ConditionItem.fix_torsion(C1.torsion(C2, C3, C4))
        theta[3] = ConditionItem.fix_torsion(C2.torsion(C3, C4, O4))
        theta[4] = ConditionItem.fix_torsion(C3.torsion(C4, O4, C1))

        p, sd_p, tm, sd_tm = ConditionItem._pseudorotation_with_sd(theta)
        return p

    def distance(self, conditional_restraint, atoms):
        vector1 = []
        vector2 = []
        for restraint in conditional_restraint:
            if restraint.name in self.restraint_names:
                vector1.append(restraint.value(atoms))
                value = getattr(self, 'atoms_%s' % restraint.type)(restraint, atoms)
                vector2.append(value)
        return self.measure(vector1, vector2)


class ConditionalRestraintItem(object):
    # use slots to minimize memory footprint
    __slots__ = (
        'type',
        'name',
        'atom_names',
        '_value',
        '_sigma',
        '_regressor',
        'value_param_name',
        'value_param_atoms',
        'value_dist',
        'sigma_dist',
    )

    def __init__(self, restraint_type, name, atom_names, value, sigma, value_dist=None, sigma_dist=None, value_function_name=None, value_param_def=None):
        self.type = restraint_type
        self.name = name
        self.atom_names = atom_names
        self._value = value
        self._sigma = sigma
        self._regressor = None
        if value_function_name:
            self._regressor = joblib.load(self.regressor_absolute_path(value_function_name))

        self.value_param_name = None
        self.value_param_atoms = None
        if value_param_def:
            self.value_param_name = value_param_def[0]
            self.value_param_atoms = value_param_def[1]
            if self.value_param_name not in ('tau_max', 'torsion_chi'):
                raise Exception('Unknown parameter for functional relationship')

        self.value_dist = value_dist
        self.sigma_dist = sigma_dist

    def regressor_absolute_path(self, value_function_filename):
        return os.path.join(os.path.abspath(os.path.dirname(__file__)), 'lib', 'regressors', value_function_filename)

    def calc_param(self, param_atoms):
        if self.value_param_name == 'tau_max':
            C1 = param_atoms[0]
            C2 = param_atoms[1]
            C3 = param_atoms[2]
            C4 = param_atoms[3]
            O4 = param_atoms[4]

            theta = [0.0, 0.0, 0.0, 0.0, 0.0]
            theta[0] = ConditionItem.fix_torsion(C4.torsion(O4, C1, C2))
            theta[1] = ConditionItem.fix_torsion(O4.torsion(C1, C2, C3))
            theta[2] = ConditionItem.fix_torsion(C1.torsion(C2, C3, C4))
            theta[3] = ConditionItem.fix_torsion(C2.torsion(C3, C4, O4))
            theta[4] = ConditionItem.fix_torsion(C3.torsion(C4, O4, C1))

            p, sd_p, tm, sd_tm = ConditionItem._pseudorotation_with_sd(theta)
            return tm

        elif self.value_param_name == 'torsion_chi':
            atom1 = param_atoms[0]
            atom2 = param_atoms[1]
            atom3 = param_atoms[2]
            atom4 = param_atoms[3]
            chi = ConditionItem.fix_torsion(atom1.torsion(atom2, atom3, atom4))
            return chi

    def map_atoms(self, atom_map, atom_names):
        return [atom_map[atom_name] for atom_name in atom_names]

    def value_sigma(self, atom_map):
        if self.value_param_atoms:
            param_atoms = self.map_atoms(atom_map, self.value_param_atoms)
            regression_param = self.calc_param(param_atoms)
            value, sigma = self._regressor.predict([[regression_param]], return_std=True)
            return value[0], sigma[0]
        return self._value, self._sigma

    def get_restraint(self, atom_map, conditional_restraint=None):
        atoms = self.map_atoms(atom_map, self.atom_names)
        # map params atoms and calculate param
        restraint_value, restraint_sigma = self.value_sigma(atom_map)
        return Restraint(self.type, atoms, restraint_value, restraint_sigma, self.value_dist, self.sigma_dist, self.name, conditional_restraint.name)

    def value(self, atom_map):
        return self.value_sigma(atom_map)[0]

    def sigma(self, atom_map):
        return self.value_sigma(atom_map)[1]


class ConditionItem(object):
    # use slots to minimize memory footprint
    __slots__ = (
        'multiplier',
        'type',
        'name',
        'atom_names',
        '_value',
        '_sigma',
    )

    multiplier = 4

    def __init__(self, condition_type, name, atom_names, value, sigma):
        if condition_type not in ('torsion', 'pseudorotation'):
            raise Exception('Unknown condition type')
        self.type = condition_type
        self.name = name
        self.atom_names = atom_names
        self._value = self.fix_torsion(value)
        self._sigma = sigma

    def value(self, atom_map=None):
        return self._value

    def sigma(self, atom_map=None):
        return self._sigma

    @classmethod
    def fix_torsion(cls, value):
        """
        Normalize torsion angle to be between [-180,180]
        :param value: torsion angle in deg
        :return:  torsion angle in deg  between [-180,180]
        """
        if value > 180:
            return value - ((int(value)//180 % 2) + int(value)//360)*360
        elif value < -180:
            return value - ((int(value)//180 % 2) + int(value)//360)*360
        return value

    @classmethod
    def _pseudorotation_with_sd(cls, theta):
        """
        Calculate pseudorotation angle
        :param theta: list of theta torsion values, for example theta[0] = torsion(C4', O4' C1', C2')
        :return: P in deg, standard deviation of P, Tm, standard deviation of Tm
        """
        # the initial definition  is Theta(1) = C1-C2-C3-C4, Theta(2) = C2-C3-C4-O4, etc.
        _theta = [theta[2], theta[3], theta[4], theta[0], theta[1]]

        sum_sin = 0.0
        sum_cos = 0.0

        for i_t, t in enumerate(_theta):
            x = 0.8 * math.pi * i_t
            sum_sin += t * math.sin(x)
            sum_cos += t * math.cos(x)

        P_deg = math.degrees(math.atan2(-sum_sin, sum_cos))

        if P_deg < 0.0:
            P_deg += 360.0

        P_rad = math.radians(P_deg)
        Tm = 0.4 * (math.cos(P_rad) * sum_cos - math.sin(P_rad) * sum_sin)

        ST = 0.0
        Thc = [0.0, 0.0, 0.0, 0.0, 0.0]

        for i_t, t in enumerate(_theta):
            Thc[i_t] = Tm * math.cos(P_rad+(0.8 * math.pi * i_t))
            d = t - Thc[i_t]
            ST += d * d

        sd_Tm = math.sqrt(0.4 * ST / 3.0)
        sd_P = sd_Tm / math.radians(Tm)
        return P_deg, sd_P, Tm, sd_Tm

    @classmethod
    def _pseudorotation(cls, theta):
        """
        Calculate pseudorotation angle based on
        http://dl.taq.ir/science/principles_of_Nucleic_Acid_Structure_Neidle.pdf page 29
        :param theta: list of theta torsion values, for example theta[0] = torsion(C4', O4' C1', C2')
        :return: P in deg, Tm
        """
        P_tan_numerator =  (theta[4]+theta[1]) - (theta[3]+theta[0])
        P_tan_denominator = 2.0 * theta[2] * (math.sin(math.radians(36.0)) + math.sin(math.radians(72.0)))
        P_rad = math.atan2(P_tan_numerator, P_tan_denominator)
        Tm = theta[2] / math.cos(P_rad)
        return math.degrees(P_rad), Tm

    def check_condition(self, atoms):
        if self.type == 'torsion':
            if len(self.atom_names) == 4:
                atom0 = atoms[self.atom_names[0]]
                atom1 = atoms[self.atom_names[1]]
                atom2 = atoms[self.atom_names[2]]
                atom3 = atoms[self.atom_names[3]]

                torsion = atom0.torsion(atom1, atom2, atom3)
                torsion = self.fix_torsion(torsion)

                low_bondary = self.value() - self.multiplier * self.sigma()
                high_bondary = self.value() + self.multiplier * self.sigma()
                for x in [0, 360, -360]:
                    if low_bondary <= torsion + x <= high_bondary:
                        #print 'condition true:', self.atom_names[0], self.atom_names[1], self.atom_names[2], self.atom_names[3], self.value-self.multiplier*self.sigma <= torsion+x <= self.value+self.multiplier*self.sigma, self.value-self.multiplier*self.sigma, torsion, self.value+self.multiplier*self.sigma
                        return True
                #print 'condition false:', self.atom_names[0], self.atom_names[1], self.atom_names[2], self.atom_names[3], self.value-self.multiplier*self.sigma <= torsion <= self.value+self.multiplier*self.sigma, self.value-self.multiplier*self.sigma, torsion, self.value+self.multiplier*self.sigma
                return False
            else:
                raise Exception('Wrong number of atoms for torsion condition')
        elif self.type == 'pseudorotation':
            if len(self.atom_names) == 5:
                if self.atom_names != ["C1'", "C2'", "C3'", "C4'", "O4'"]:
                    raise Exception('Wrong order of atoms for pseudorotation condition')

                C1 = atoms[self.atom_names[0]]
                C2 = atoms[self.atom_names[1]]
                C3 = atoms[self.atom_names[2]]
                C4 = atoms[self.atom_names[3]]
                O4 = atoms[self.atom_names[4]]

                theta = [0.0, 0.0, 0.0, 0.0, 0.0]
                theta[0] = self.fix_torsion(C4.torsion(O4, C1, C2))
                theta[1] = self.fix_torsion(O4.torsion(C1, C2, C3))
                theta[2] = self.fix_torsion(C1.torsion(C2, C3, C4))
                theta[3] = self.fix_torsion(C2.torsion(C3, C4, O4))
                theta[4] = self.fix_torsion(C3.torsion(C4, O4, C1))

                p, sd_p, tm, sd_tm = self._pseudorotation_with_sd(theta)

                #print("{} {} {} {} {}, p={}({}) tm={}({})".format(
                #    self.atom_names[0], self.atom_names[1], self.atom_names[2], self.atom_names[3], self.atom_names[4],
                #    p, sd_p, tm, sd_tm)
                #)

                # multiplier is 4 so possible source of error
                # TODO
                low_bondary = self.value() - self.multiplier * self.sigma()
                high_bondary = self.value() + self.multiplier * self.sigma()
                for x in [0, 360, -360]:
                    if low_bondary <= p + x <= high_bondary:
                        return True
                return False
            else:
                raise Exception('Wrong number of atoms for pseudorotation condition')
        raise Exception('Unknown condition type')


class ConditionalRestraint(object):
    # use slots to minimize memory footprint
    __slots__ = (
        'name',
        'conditions',
        'restraints',
    )

    def __init__(self, name, conditions, restraints):
        create_condition = lambda condition: condition if isinstance(condition, ConditionItem) else ConditionItem(*condition)
        create_restraint = lambda restraint: restraint if isinstance(restraint, ConditionalRestraintItem) else ConditionalRestraintItem(*restraint)

        self.name = name
        self.conditions = [create_condition(con) for con in conditions]
        self.restraints = [create_restraint(res) for res in restraints]

    def is_default(self):
        if len(self.conditions) == 0 or self.name == 'default':
            return True
        return False

    def check_conditions(self, atoms):
        if self.is_default():
            return True
        for condition in self.conditions:
            if condition.check_condition(atoms) is False:
                return False
        return True

    def get_restraints(self, atoms):
        printable_restraints = []
        for restraint in self.restraints:
            printable_restraints.append(restraint.get_restraint(atoms, self))
        return printable_restraints


class ConditionalRestraintList(list):
    def __init__(self, data=()):
        super(ConditionalRestraintList, self).__init__()
        create_conditional_restraint = lambda obj: obj if isinstance(obj, ConditionalRestraint) else \
            ConditionalRestraint(**obj)
        for item in data:
            self.append(create_conditional_restraint(item))

    def __getslice__(self, i, j):
        return ConditionalRestraintList(list.__getslice__(self, i, j))

    def get_feasible(self, atoms):
        feasible = ConditionalRestraintList([])
        for conditional_restraint in iter(self):
            if conditional_restraint.check_conditions(atoms) is True:
                feasible.append(conditional_restraint)
        return feasible

    def any_default(self):
        return any([conditional_restraint.is_default() for conditional_restraint in self])

    def remove_default(self):
        i_to_delete = []
        for i, conditional_restraint in enumerate(self):
            if conditional_restraint.is_default():
                i_to_delete.append(i)
        for i in reversed(i_to_delete):
            self.pop(i)

    def get_default(self):
        default = None
        for i, conditional_restraint in enumerate(self):
            if conditional_restraint.is_default():
                return conditional_restraint
        return default

    def find_closest(self, atoms, distance_measure, variable):

        if len(self) == 0:
            return None

        if self.any_default() and len(self) > 1:
            self.remove_default()

        min_distance = 999999999.0
        min_distance_i = None
        for i, conditional_restraint in enumerate(self):
            distance = distance_measure.distance(getattr(conditional_restraint, variable), atoms)
            #if "PO4==AS" in conditional_restraint.name:
            #    print(i, conditional_restraint.name, distance, [str(atom) for atom in atoms.values()])
            if distance < min_distance:
                min_distance = distance
                min_distance_i = i

        if min_distance_i is not None:
            return self[min_distance_i]

        return None

    def find_restraint_closest(self, atoms, distance_measure):
        return self.find_closest(atoms, distance_measure, 'restraints')

    def find_condition_closest(self, atoms, distance_measure):
        return self.find_closest(atoms, distance_measure, 'conditions')


class Restraint(object):
    # use slots to minimize memory footprint
    __slots__ = (
        'type',
        'atoms',
        'value',
        'sigma',
        'value_dist',
        'sigma_dist',
        'name',
        'condition_name',
    )

    def __init__(self, restraint_type, atoms, value, sigma, value_dist=None, sigma_dist=None, name=None, condition_name=None):
        self.type = restraint_type
        self.atoms = atoms
        self.value = value
        if isinstance(self.value, str):
            raise Exception('Value should be float type')
        self.sigma = sigma
        if isinstance(self.sigma, str):
            raise Exception('Value should be float type')
        self.value_dist = value_dist
        self.sigma_dist = sigma_dist
        self.name = name
        self.condition_name = condition_name


class AtomGroupCache:
    # use slots to minimize memory footprint
    __slots__ = (
        '_atoms',
        '_neighbours',
    )

    def __init__(self):
        self._atoms = []
        self._neighbours = []

    def add_atom(self, atom_idx):
        self._atoms.append(atom_idx)

    def iter_atoms(self):
        return self._atoms

    def find_atom(self, chain_id, atom_name, res_id):
        for atom in self._atoms:
            if atom.atom_name == atom_name and atom.res_id == res_id and atom.chain_id == chain_id:
                return atom
        return None

    def add_neighbour(self, neighbour):
        self._neighbours.append(neighbour)

    def iter_neighbour(self):
        return self._neighbours

    def find_neighbour(self, chain_id, atom_name, res_id, alt_loc):
        for atom in self._neighbours:
            if (
                atom.atom_name == atom_name and
                atom.res_id == res_id and
                atom.chain_id == chain_id and
                (atom.alt_loc==alt_loc or (atom.alt_loc=='' and alt_loc!='') or (alt_loc=='' and atom.alt_loc!=''))
            ):
                return atom
        return None


class MonomerRestraintGroup(object):
    # use slots to minimize memory footprint
    __slots__ = (
        'name',
        'res_names',
        'atom_labels',
        '_valid_atom_labels',
        'res_numbers',
        'conditions',
        'disallowed_conditions',
        'distance_measure',
        'condition_measure',
        'restraints',
        'atoms',
        'groups',
        '_registered_res_id'
    )

    def __init__(
            self,
            name,
            res_names,
            atom_labels,
            res_numbers,
            conditions,
            disallowed_conditions,
            restraints,
            distance_measure,
            condition_measure
    ):
        self.name = name
        self.res_names = res_names
        self.atom_labels = atom_labels
        self._valid_atom_labels = atom_labels.keys()
        self.res_numbers = res_numbers
        self.conditions = conditions
        self.disallowed_conditions = disallowed_conditions
        self.distance_measure = DistanceMeasure(**distance_measure)
        self.condition_measure = DistanceMeasure(**condition_measure)

        if isinstance(restraints, ConditionalRestraintList):
            self.restraints = restraints
        else:
            self.restraints = ConditionalRestraintList(restraints)
        self.atoms = []
        # dict (chain_id, res_id-modifier, altloc): [] - lists with indexes to atoms
        self.groups = defaultdict(AtomGroupCache)

        # keep the res if that matches the pdb_code
        self._registered_res_id = defaultdict(set)

    def register_valid_res_id(self, chain_id, res_id):
        self._registered_res_id[chain_id.strip()].add(int(res_id))

    def is_registered_res_id_or_neighbour(self, chain_id, res_id):
        _res_id = int(res_id)
        registered_res_id = self._registered_res_id.get(chain_id.strip(), set())
        return (_res_id in registered_res_id) or (_res_id+1 in registered_res_id) or (_res_id-1 in registered_res_id)

    def is_valid_res_name(self, res_name):
        return res_name.strip() in self.res_names

    def is_valid_atom_name(self, atom_name):
        return atom_name.strip() in self._valid_atom_labels

    def any_atom_in_disallowed_atoms(self, group, chain_id, res_id, alt_id):
        for atom_name_1, atom_name_2, dist, res_id_mod_1, res_id_mod_2 in self.disallowed_conditions:
            # self.atoms can be inefficient
            # we should get it from the self.group[group_key]
            # to do so we need to refactor self.group dict into proper class
            atom_1 = group.find_neighbour(chain_id, atom_name_1, res_id + res_id_mod_1, alt_id)
            atom_2 = group.find_neighbour(chain_id, atom_name_2, res_id + res_id_mod_2, alt_id)
            if atom_1 is not None and atom_2 is not None and atom_1.dist(atom_2) <= dist:
                return True
            #print("Disallowed atoms not found", self.name, chain_id, atom_name_1, res_id + res_id_mod_1, alt_id, chain_id, atom_name_2, res_id + res_id_mod_2, alt_id)
        return False

    def add_atom(self, chain_id, res_id, res_name, atom_name, alt_loc, atom_xyz):
        if self.is_registered_res_id_or_neighbour(chain_id, res_id):
            atom = Atom(chain_id, res_id, res_name, atom_name, alt_loc, atom_xyz)
            self.atoms.append(atom)

    def create_res_groups(self):
        preliminary_groups = defaultdict(AtomGroupCache)

        for chain_id, registered_res_id in self._registered_res_id.iteritems():
            for res_id in registered_res_id:
                for i_atom, atom in enumerate(self.atoms):
                    key = "{}_{:05d}".format(chain_id, res_id)
                    if (
                        chain_id == atom.chain_id and
                        atom.atom_name in self.res_numbers and
                        res_id == atom.res_id - self.res_numbers[atom.atom_name]
                    ):
                        preliminary_groups[key].add_atom(atom)
                    if (
                        chain_id == atom.chain_id and
                        abs(res_id - atom.res_id) <= 1
                    ):
                        preliminary_groups[key].add_neighbour(atom)

        for key, preliminary_group in preliminary_groups.iteritems():
            p_group_atoms = preliminary_group.iter_atoms()
            locs = set([atom.alt_loc for atom in p_group_atoms])
            locs.discard('')

            if len(locs) > 0:
                # at least one atom with alt_loc code
                for loc in locs:
                    #create new key
                    key_alt = "{}_{}".format(key, loc)
                    for atom in p_group_atoms:
                        atom_alt_loc = atom.alt_loc
                        if atom_alt_loc in (loc, ''):
                            self.groups[key_alt].add_atom(atom)
                    self.groups[key_alt]._neighbours = preliminary_group._neighbours
            else:
                # no alt_loc codes
                self.groups[key] = preliminary_group

    def _print_groups(self):
        print "# PRINT GROUPS {}".format(self.name)
        for key, group in self.groups.iteritems():
            for atom in group.iter_atoms():
                print self.name, key, atom.chain_id, atom.res_id, atom.res_name, atom.atom_name, atom.alt_loc, atom.atom_xyz
        print "# PRINT GROUPS finished"

    def is_valid_atom_group(self, group_key, group):
        split_key = group_key.split("_")
        chain_id = split_key[0].strip()
        res_id = int(split_key[1])
        alt_id = '' if len(split_key) < 3 else split_key[2].strip()

        if self.any_atom_in_disallowed_atoms(group, chain_id, res_id, alt_id):
            return False

        for atom_name_1, atom_name_2, dist, res_id_mod_1, res_id_mod_2 in self.conditions:
            # if it is atom1.alt_id == '' and atom2.alt_id=='A', atom2.alt_id=='B' there is a problem
            atom_1 = group.find_neighbour(chain_id, atom_name_1, res_id + res_id_mod_1, alt_id)
            atom_2 = group.find_neighbour(chain_id, atom_name_2, res_id + res_id_mod_2, alt_id)

            if atom_1 is None or atom_2 is None or atom_1.dist(atom_2) > dist:
                print "{} {} {}: non valid dist: {}, atom1: {} {} {}, atom2: {} {} {}".format(
                    self.name,
                    chain_id,
                    res_id,
                    'at least one missing' if atom_1 is None or atom_2 is None else atom_1.dist(atom_2),
                    atom_1.chain_id if atom_1 else '',
                    atom_1.res_id if atom_1 else '',
                    atom_name_1,
                    atom_2.chain_id if atom_2 else '',
                    atom_2.res_id if atom_2 else '',
                    atom_name_2,
                )
                return False
        return True

    def validate_links(self, verbose=0):
        """
        Deletes groups that are not bonded correctly
        """
        if verbose > 0:
            print '#'*60
            print "# SEARCHING for {} group/monomer".format(self.name)
        for group_key, group in sorted(self.groups.iteritems()):
            if not self.is_valid_atom_group(group_key, group):
                del self.groups[group_key]
                if verbose > 0:
                    print "#     {} ignoring".format(group_key)
            else:
                if verbose > 0:
                    print "#     {} recognized as {}".format(group_key, self.name)
        if verbose > 0:
            print "# SEARCHING finished".format()

    def atom_restraints(self, verbose=0):
        result_restraints = []
        if verbose > 0:
            print '# ANALYZING {} group/monomer'.format(self.name)
        for atom_group_key, atom_group in sorted(self.groups.iteritems()):
            atoms = {atom.atom_name: atom for atom in atom_group.iter_atoms()}
            if verbose > 1:
                print "#     {} atoms for {} {}".format(atom_group_key, self.name, list(atoms.keys()))

            feasible_restraints = self.restraints.get_feasible(atoms)
            if verbose > 1:
                print "#     {} feasible_restraints for {} {}".format(atom_group_key, self.name, [
                    fc.name for fc in feasible_restraints
                ])

            closest_restraint = feasible_restraints.find_restraint_closest(atoms, self.distance_measure)

            if closest_restraint is not None:
                if verbose > 0:
                    print "#     {} recognized as {} {}".format(atom_group_key, self.name, closest_restraint.name)
            else:
                if verbose > 0:
                    print "#     {} not feasible".format(atom_group_key)

                closest_restraint = self.restraints.get_default()

                if closest_restraint is None:
                    closest_restraint = self.restraints.find_condition_closest(atoms, self.condition_measure)
                    if verbose > 0:
                        print "#             closest to {} {}".format(self.name, closest_restraint.name)
                else:
                    if verbose > 0:
                        print "#             default used"
                        print "#             set default {} {}".format(self.name, closest_restraint.name)

            result_restraints.extend(closest_restraint.get_restraints(atoms))
        if verbose > 0:
            print '# ANALYZING finished'
            print '#'*60
        return result_restraints

    def prepare_restraints(self, verbose=0):
        self.create_res_groups()
        if verbose > 0:
            self._print_groups()
        self.validate_links(verbose=verbose)

        restraints = self.atom_restraints(verbose=verbose)
        return restraints


def save(stream, restraint_groups, restraint_text_all, printer_cls, in_filename):
    printer_cls.save_info(stream, restraint_groups, version=VERSION)
    printer_cls.save_input_filename(stream, in_filename)

    if len(restraint_text_all) == 0:
        print >> stream, "# There were no restraints to be created based on the submitted PDB file"
    else:
        header = printer_cls.header()
        if header != '':
            print >> stream, header

        print >> stream, "\n".join(restraint_text_all)

        footer = printer_cls.footer()
        if footer != '':
            print >> stream, footer


def print_info_about_restraints(all_restraints):
    res_id_name_map = defaultdict(set)
    atom_cache = defaultdict(set)
    restraint_cache = defaultdict(set)

    for restraint in all_restraints:
        for atom in restraint.atoms:
            key = "%s_%05d_%s%s" % (atom.chain_id, atom.res_id, atom.res_name, ("_%s" % atom.alt_loc if atom.alt_loc else ''))
            value = restraint.condition_name
            res_id_name_map[key].add(value)

            cache_key = "%s_%s" % (key, value)
            atom_cache[cache_key].add(atom.atom_name+("^%s" % atom.alt_loc if atom.alt_loc else ''))
            restraint_cache[cache_key].add(restraint.name)

    for key in sorted(res_id_name_map.keys()):
        print(key)
        for val in sorted(res_id_name_map[key]):
            cache_key = "%s_%s" % (key, val)
            atoms = atom_cache[cache_key]
            rest = restraint_cache[cache_key]
            print("\t%s" % (val))
            print("\t\tAtoms: %s" % (sorted(atoms)))
            print("\t\tDist: %s" % (sorted((res for res in rest if res.startswith('d')))))
            print("\t\tAngles: %s" % (sorted((res for res in rest if res.startswith('a')))))


def parse_pdb(in_pdb, restraint_groups, allowed_restraint_groups, out_filename, printer_cls, lines=False, source_info=None):
    if lines is False and (in_pdb.endswith(".res") or in_pdb.endswith(".ins")):
        raise Exception("Shelx format files not supported")
    else:
        in_filename = source_info if lines else in_pdb
        data_pdb = iotbx.pdb.input(lines=flex.split_lines(in_pdb), source_info=source_info) if lines else iotbx.pdb.input(file_name=in_pdb)
        pdb_hierarchy = data_pdb.construct_hierarchy()

        printer_cls.validate(pdb_hierarchy)

        for model in pdb_hierarchy.models():
            for chain in model.chains():
                for residue_group in chain.residue_groups():
                    for atom_group in residue_group.atom_groups():
                        # if necessary only for speed optimization
                        for restraint in restraint_groups:
                            if restraint.is_valid_res_name(atom_group.resname):
                                restraint.register_valid_res_id(chain.id, residue_group.resid())


        for model in pdb_hierarchy.models():
            for chain in model.chains():
                for residue_group in chain.residue_groups():
                    for atom_group in residue_group.atom_groups():
                        # if necessary only for speed optimization
                        for restraint in restraint_groups:
                            if restraint.is_registered_res_id_or_neighbour(chain.id, residue_group.resid()):
                                altloc = atom_group.altloc.strip()
                                for atom in atom_group.atoms():
                                    restraint.add_atom(
                                        chain.id, residue_group.resid(), atom_group.resname, atom.name, altloc, atom.xyz
                                    )

    all_restraints = []
    for restraint_group in restraint_groups:
        all_restraints.extend(restraint_group.prepare_restraints())

    # print_info_about_restraints(all_restraints)

    restraint_text_all = []
    restraint_text = printer_cls.print_restraints(all_restraints, allowed_restraint_groups)
    if len(restraint_text) > 0:
        restraint_text_all.append(restraint_text)

    if type(out_filename) == str or type(out_filename) == unicode:
        with open(out_filename, 'w') as res_file:
            save(res_file, allowed_restraint_groups, restraint_text_all, printer_cls, in_filename)
    else:
        save(out_filename, allowed_restraint_groups, restraint_text_all, printer_cls, in_filename)


def create_monomer_group(module, prefix, name):
    return MonomerRestraintGroup(
        name,
        getattr(module, '%s_PDB_CODES' % prefix),
        getattr(module, '%s_ATOM_NAMES' % prefix),
        getattr(module, '%s_ATOM_RES' % prefix),
        getattr(module, '%s_REQUIRED_CONDITION' % prefix),
        getattr(module, '%s_DISALLOWED_CONDITION' % prefix, {}),
        getattr(module, '%s_RESTRAINTS' % prefix),
        getattr(module, '%s_DISTANCE_MEASURE' % prefix),
        getattr(module, '%s_CONDITION_DISTANCE_MEASURE' % prefix),
    )


def load_restraints_lib(po4=True, po4terminal=True, bases=True, isobases=True, ribose_deoxyribose=True,
                        ribose_deoxyribose_terminal=True):
    #TODO
    #should dynamically load from lib
    restraint_list = []
    if po4:
        restraint_list.append(create_monomer_group(PO4, 'PO4', 'PO4'))
    if po4terminal:
        restraint_list.append(create_monomer_group(PO4_terminal_C5, 'PO4_5_TERMINAL', 'PO4_terminal_C5'))
        restraint_list.append(create_monomer_group(PO4_terminal_C3, 'PO4_3_TERMINAL', 'PO4_terminal_C3'))
    if bases:
        for prefix in ('ADENINE', 'GUANINE', 'URACIL', 'THYMINE', 'CYTOSINE'):
            restraint_list.append(create_monomer_group(nucleic_acid_bases, prefix, prefix))
    if isobases:
        for prefix in ('ISOCYTOSINE', 'ISOGUANINE'):
            restraint_list.append(create_monomer_group(nucleic_acid_isobases, prefix, prefix))
    if ribose_deoxyribose:
        modules = [deoxyribose_purine, deoxyribose_pyrimidine, ribose_purine, ribose_pyrimidine]
        prefixes = ['DEOXYRIBOSE_PURINE', 'DEOXYRIBOSE_PYRIMIDINE', 'RIBOSE_PURINE', 'RIBOSE_PYRIMIDINE']
        for module, prefix in zip(modules, prefixes):
            for group in [
                'CHI_CONFORMATION',
                'GAMMA',
                'CONFORMATION',
                'BASE_FUNC_OF_TORSION_CHI',
                'ALL_FUNC_OF_TORSION_CHI',
                'SUGAR_CONFORMATION_FUNC_OF_TAU_MAX',
                'CHI',
                'CHI_GAMMA',
                'ALL',
                'SUGAR'
            ]:
                prefix_group = '{}_{}'.format(prefix, group)
                restraint_list.append(create_monomer_group(module, prefix_group, prefix_group))
    if ribose_deoxyribose_terminal:
        modules = [deoxyribose_purine_terminal_C3, deoxyribose_pyrimidine_terminal_C3,
                   ribose_purine_terminal_C3, ribose_pyrimidine_terminal_C3,
                   ]
        prefixes = ['DEOXYRIBOSE_PURINE_TERMINAL_C3', 'DEOXYRIBOSE_PYRIMIDINE_TERMINAL_C3',
                    'RIBOSE_PURINE_TERMINAL_C3', 'RIBOSE_PYRIMIDINE_TERMINAL_C3',
                    ]
        for module, prefix in zip(modules, prefixes):
            for group in [
                'CHI_CONFORMATION',
                'GAMMA',
                'CONFORMATION',
                'BASE_FUNC_OF_TORSION_CHI',
                'ALL_FUNC_OF_TORSION_CHI',
                'SUGAR_CONFORMATION_FUNC_OF_TAU_MAX',
                'CHI',
                'CHI_GAMMA',
                'ALL',
                'SUGAR',
                'SUGAR_CONFORMATION'
            ]:
                prefix_group = '{}_{}'.format(prefix, group)
                restraint_list.append(create_monomer_group(module, prefix_group, prefix_group))

        modules = [deoxyribose_purine_terminal_C5, deoxyribose_pyrimidine_terminal_C5,
                   ribose_purine_terminal_C5, ribose_pyrimidine_terminal_C5,
                   ]
        prefixes = ['DEOXYRIBOSE_PURINE_TERMINAL_C5', 'DEOXYRIBOSE_PYRIMIDINE_TERMINAL_C5',
                    'RIBOSE_PURINE_TERMINAL_C5', 'RIBOSE_PYRIMIDINE_TERMINAL_C5',
                    ]
        for module, prefix in zip(modules, prefixes):
            for group in [
                'CHI_CONFORMATION',
                'GAMMA',
                'CONFORMATION',
                'BASE_FUNC_OF_TORSION_CHI',
                'ALL_FUNC_OF_TORSION_CHI',
                'SUGAR_CONFORMATION_FUNC_OF_TAU_MAX',
                'CHI',
                'CHI_GAMMA',
                'ALL',
                'SUGAR'
            ]:
                prefix_group = '{}_{}'.format(prefix, group)
                restraint_list.append(create_monomer_group(module, prefix_group, prefix_group))

    return restraint_list


def run():
    if len(sys.argv) > 3:
        printer = sys.argv[1]
        in_pdb = sys.argv[2]
        out_filename = sys.argv[3]
    elif len(sys.argv) > 2:
        printer = 'Refmac'
        in_pdb = sys.argv[1]
        out_filename = sys.argv[2]
    else:
        print 'usage: restraints.py [printer] in.pdb restraints.txt'
        print '       Printer is one of Refmac, Phenix, Shelxl, Csv. Default=Refmac'
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
        print "Unknown printer {}, should be one of Refmac, Phenix, Shelxl".format(printer)
        return

    restraint_list = load_restraints_lib()
    parse_pdb(in_pdb, restraint_list, restraint_list, out_filename, printer_cls)


if __name__ == "__main__":
    run()
