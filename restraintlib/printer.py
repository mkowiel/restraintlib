# coding: utf-8
from __future__ import print_function

import logging
import math
import os
from datetime import datetime

__author__ = "Marcin Kowiel, Dariusz Brzezinski"


class PrinterBase(object):

    def __init__(self, override_sigma=True):
        self.override_sigma = override_sigma

    @classmethod
    def validate(cls, hierarchy):
        return None

    @classmethod
    def header(cls):
        return ''

    @classmethod
    def footer(cls):
        return ''

    def angle_sigma_value(self, value):
        return value

    def get_angle(self, restraint, all_restraints):
        return ''

    def dist_sigma_value(self, value):
        return value

    def get_dist(self, restraint, all_restraints):
        return ''

    def print_restraints(self, restraints, allowed_restraint_groups):
        lines = []

        allowed_condition_name = [
            conditional_restraint.name
            for group in allowed_restraint_groups
            for conditional_restraint in group.restraints
        ]

        for restraint in restraints:
            # for shelx printer we need dist to print angles, they can be defined in other lib
            if restraint.condition_name in allowed_condition_name:
                line = ''
                if restraint.type == 'angle':
                    line = self.get_angle(restraint, restraints)
                elif restraint.type == 'dist':
                    line = self.get_dist(restraint, restraints)

                # usuwanie duplikujacych sie wynikow dla alt_loc (blank-blank) w przypadku gdy mamy A i B
                if line not in lines:
                    lines.append(line)

        return "\n".join(lines)

    @classmethod
    def header_comment_begin(cls):
        return '#'

    @classmethod
    def header_comment_end(cls):
        return '#'

    @classmethod
    def save_info(cls, stream, restraint_list, version):
        libs_lines = [""]
        for res in restraint_list:
            appended_line = "{}{}, ".format(libs_lines[-1], res.name)
            if len(appended_line) < 70:
                libs_lines[-1] = appended_line
            else:
                libs_lines.append("{}, ".format(res.name))

        comm_begin = cls.header_comment_begin()
        comm_end = cls.header_comment_end()

        print("{}#########################################################################{}".format(comm_begin, comm_end), file=stream)
        print("{}                   RestraintLib version {}                        {}".format(comm_begin, version, comm_end), file=stream)
        print("{}#########################################################################{}".format(comm_begin, comm_end), file=stream)
        print("{}                                                                         {}".format(comm_begin, comm_end), file=stream)
        print("{} M.Kowiel, D.Brzezinski, M.Jaskolski (2016)                              {}".format(comm_begin, comm_end), file=stream)
        print("{} Conformation-dependent restraints for polynucleotides:                  {}".format(comm_begin, comm_end), file=stream)
        print("{} I. Clustering of the geometry of the phosphodiester group.              {}".format(comm_begin, comm_end), file=stream)
        print("{} Nucleic Acids Res. 44, 8479–8489.                                       {}".format(comm_begin, comm_end), file=stream)
        print("{} https://doi.org/10.1093/nar/gkw717                                      {}".format(comm_begin, comm_end), file=stream)
        print("{}                                                                         {}".format(comm_begin, comm_end), file=stream)
        print("{} M.Gilski, J.Zhao, M.Kowiel, D.Brzezinski, D.H.Turner, M.Jaskolski (2019){}".format(comm_begin, comm_end), file=stream)
        print("{} Accurate geometrical restraints for Watson-Crick base pairs.            {}".format(comm_begin, comm_end), file=stream)
        print("{} Acta Cryst. B75, 235-245.                                               {}".format(comm_begin, comm_end), file=stream)
        print("{} https://doi.org/10.1107/S2052520619002002                               {}".format(comm_begin, comm_end), file=stream)
        print("{}                                                                         {}".format(comm_begin, comm_end), file=stream)
        print("{} M.Kowiel, D.Brzezinski, M. Gilski, M.Jaskolski (2020)                   {}".format(comm_begin, comm_end), file=stream)
        print("{} Conformation-dependent restraints for polynucleotides: The sugar moiety.{}".format(comm_begin, comm_end), file=stream)
        print("{} Nucleic Acids Res. 48, 962–973.                                         {}".format(comm_begin, comm_end), file=stream)
        print("{} https://doi.org/10.1093/nar/gkz1122                                     {}".format(comm_begin, comm_end), file=stream)
        print("{}                                                                         {}".format(comm_begin, comm_end), file=stream)
        print("{}#########################################################################{}".format(comm_begin, comm_end), file=stream)
        print("{}                                                                         {}".format(comm_begin, comm_end), file=stream)
        print("{} Restraints for:                                                         {}".format(comm_begin, comm_end), file=stream)
        for line in libs_lines:
            print("{}  {}{} {}".format(comm_begin, line, ' ' * (70 - len(line)), comm_end), file=stream)
        print("{}                                                                         {}".format(comm_begin, comm_end), file=stream)
        print("{}#########################################################################{}".format(comm_begin, comm_end), file=stream)

    @classmethod
    def save_input_filename(cls, stream, in_filename):
        comm_begin = cls.header_comment_begin()
        comm_end = cls.header_comment_end()

        print("{} Time of creation: {}Z".format(comm_begin, datetime.utcnow().isoformat()), file=stream)
        print("{} Input filename: {}".format(comm_begin, os.path.basename(in_filename)), file=stream)
        print("{}#########################################################################{}".format(comm_begin, comm_end), file=stream)


class ShelxPrinter(PrinterBase):
    @classmethod
    def header_comment_begin(cls):
        return 'REM '

    @classmethod
    def header_comment_end(cls):
        return ''

    @classmethod
    def _get_alt_loc(cls, atom):
        return '^{}'.format(atom.alt_loc).lower().strip() if str(atom.alt_loc).strip() != '' else ''

    @classmethod
    def _get_resid(cls, atom):
        resid = str(atom.res_id).strip()
        chainid = str(atom.chain_id).strip()
        return resid if chainid == '' else ("%s:%s" % (chainid, resid))

    @classmethod
    def _comment(cls, restraint):
        if restraint.condition_name is not None and restraint.name is not None:
            if len(restraint.condition_name) < 64:
                restraint_condition_name_line = 'REM Restraint {}'.format(restraint.condition_name)
            else:
                restraint_condition_name_line = 'REM Restraint\nREM {}'.format(restraint.condition_name)
            if restraint.type == 'angle':
                return '{}\nREM {} {} {:.1f} {:.1f}'.format(
                    restraint_condition_name_line, restraint.type, restraint.name, restraint.value, restraint.sigma)
            return '{}\nREM {} {} {:.3f} {:.3f}'.format(
                restraint_condition_name_line, restraint.type, restraint.name, restraint.value, restraint.sigma)
        return ''

    @classmethod
    def _get_resid_altcode_rem(cls, atom):
        resid = cls._get_resid(atom)
        alt_code = cls._get_alt_loc(atom)
        chainid = str(atom.chain_id).strip()
        if alt_code != '' and chainid != '':
            rem = 'REM Chain {} and resid {} and name {} and altloc {} renamed to {}_{}{}'
            rem = rem.format(
                atom.chain_id,
                atom.res_id,
                atom.atom_name,
                atom.alt_loc,
                atom.atom_name,
                resid,
                alt_code,
            )
            rem.replace('  ', ' ')
        else:
            rem = ''
        return resid, alt_code, rem

    def get_dist(self, restraint, all_restraints):
        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]

        resid_0, alt_code_0, rem_0 = self._get_resid_altcode_rem(atom0)
        resid_1, alt_code_1, rem_1 = self._get_resid_altcode_rem(atom1)
        comment = self._comment(restraint)

        lines = []
        for rem in (rem_0, rem_1, comment):
            if rem != '':
                for line in rem.split('\n'):
                    lines.append(line)

        lines.append('DFIX {:.3f} {:.3f} {}_{}{} {}_{}{}'.format(
            restraint.value,
            restraint.sigma,
            atom0.atom_name,
            resid_0,
            alt_code_0,
            atom1.atom_name,
            resid_1,
            alt_code_1,
        ))
        return "\n".join(lines)

    @classmethod
    def _deg_to_dist(cls, angle_in_deg, a, b):
        return math.sqrt(a*a+b*b-2*a*b*math.cos(math.radians(angle_in_deg)))

    @classmethod
    def _find_dist(cls, atom0, atom1, all_restraints):
        # TODO
        # check is only 1 option
        for res in all_restraints:
            if res.type == 'dist':
                if (atom0 in res.atoms and atom1 in res.atoms):
                    return res.value, res.sigma
        return None, None

    @classmethod
    def _val_deg_to_dist(cls, restraint, all_restraints):
        if restraint.value_dist is not None:
            return restraint.value_dist

        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]
        atom2 = restraint.atoms[2]

        angle = restraint.value
        if angle is None:
            raise Exception("Angle should not be None")

        a, a_sig = cls._find_dist(atom0, atom1, all_restraints)
        b, b_sig = cls._find_dist(atom1, atom2, all_restraints)

        if a is None:
            for res in all_restraints:
                if res.type == 'dist':
                    res_atom0 = res.atoms[0]
                    res_atom1 = res.atoms[1]
                    if (
                        res_atom0.chain_id == atom0.chain_id and
                        res_atom0.res_id == atom0.res_id and (
                            atom0.atom_name == res_atom0.atom_name or atom0.atom_name == res_atom1.atom_name or
                            atom1.atom_name == res_atom0.atom_name or atom1.atom_name == res_atom1.atom_name
                        )
                    ):
                        print(res.condition_name, res.name, [str(atom) for atom in res.atoms])
            raise Exception("For angle restraint the relevant dist restraint should be added (atoms: {}, {})".format(
                str(atom0), str(atom1)
            ))
        elif b is None:
            for res in all_restraints:
                if res.type == 'dist':
                    res_atom0 = res.atoms[0]
                    res_atom1 = res.atoms[1]
                    if (
                        res_atom0.chain_id == atom1.chain_id and
                        res_atom0.res_id == atom1.res_id and (
                            atom1.atom_name == res_atom0.atom_name or atom1.atom_name == res_atom1.atom_name or
                            atom2.atom_name == res_atom0.atom_name or atom2.atom_name == res_atom1.atom_name
                        )
                    ):
                        print(res.condition_name, res.name, [str(atom) for atom in res.atoms])
            raise Exception("For angle restraint the relevant dist restraint should be added (atoms: {}, {})".format(
                str(atom1), str(atom2)
            ))

        return cls._deg_to_dist(angle, a, b)

    @classmethod
    def _sig_deg_to_dist(cls, restraint, all_restraints):
        if restraint.sigma_dist is not None:
            return restraint.sigma_dist

        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]
        atom2 = restraint.atoms[2]

        angle = restraint.value
        angle_sigma = restraint.sigma
        a, a_sig = cls._find_dist(atom0, atom1, all_restraints)
        b, b_sig = cls._find_dist(atom1, atom2, all_restraints)

        dist13_plus = cls._deg_to_dist(angle+0.5*angle_sigma, a, b)
        dist13_minus = cls._deg_to_dist(angle-0.5*angle_sigma, a, b)

        return dist13_plus-dist13_minus

    def get_angle(self, restraint, all_restraints):
        atom0 = restraint.atoms[0]
        atom2 = restraint.atoms[2]

        resid_0, alt_code_0, rem_0 = self._get_resid_altcode_rem(atom0)
        resid_2, alt_code_2, rem_2 = self._get_resid_altcode_rem(atom2)
        comment = self._comment(restraint)

        lines = []
        for rem in (rem_0, rem_2, comment):
            if rem != '':
                for line in rem.split('\n'):
                    lines.append(line)

        lines.append('DANG {:.3f} {:.3f} {}_{}{} {}_{}{}'.format(
            self._val_deg_to_dist(restraint, all_restraints),
            self._sig_deg_to_dist(restraint, all_restraints),
            atom0.atom_name,
            resid_0,
            alt_code_0,
            atom2.atom_name,
            resid_2,
            alt_code_2
        ))
        return "\n".join(lines)


class PhenixPrinter(PrinterBase):

    @classmethod
    def header(cls):
        return 'refinement.geometry_restraints.edits {'

    @classmethod
    def footer(cls):
        return '}'

    @classmethod
    def _get_alt_loc(cls, atom):
        return ' and altloc {}'.format(atom.alt_loc) if atom.alt_loc != '' else ''

    @classmethod
    def _atom_sel(cls, atom):
        return "chain {} and resid {} and name {}{}".format(
            atom.chain_id,
            atom.res_id,
            atom.atom_name,
            cls._get_alt_loc(atom),
        )

    def angle_sigma_value(self, value):
        """ Default sigma for angles in PHENIX is 3.0"""
        return 3.0 if self.override_sigma else value

    def dist_sigma_value(self, value):
        """ Default sigma for bond distances in PHENIX is 0.020"""
        return 0.020 if self.override_sigma else value

    def action_type(self, kind, restraint):
        return 'change'

    def atoms_with_fixed_atom_order(self, restraint):
        """
        Older Phenix requires atoms to be in the right order for action=*change
        """
        if restraint.name in ("dO1P4", "dO2P4",  "dO3P4", "dO5P4"):
            if restraint.atoms[0].atom_name != 'P':
                return reversed(restraint.atoms)
        elif restraint.name in ("dO3C3", "dO5C5"):
            if restraint.atoms[0].atom_name != 'C':
                return reversed(restraint.atoms)
        return restraint.atoms

    def _prepare_restraint(self, kind, digit_after_dot, restraint):
        lines = []
        if restraint.condition_name is not None and restraint.name is not None:
            lines.append(
                '  # Restraint {{}} {{}} {{}} {{:.{}f}} {{:.{}f}}'.format(
                    digit_after_dot, digit_after_dot).format(
                    restraint.condition_name, restraint.type, restraint.name, restraint.value, restraint.sigma
                )
            )
        lines.append("  {} {{".format(kind))
        action = self.action_type(kind, restraint)
        lines.append("    action = *{}".format(action))
        for i_atom, atom in enumerate(self.atoms_with_fixed_atom_order(restraint), start=1):
            lines.append("    atom_selection_{} = {}".format(i_atom, self._atom_sel(atom)))
        #lines.append("    symmetry_operation = None")
        keyword_ideal = 'distance_ideal' if kind == 'bond' else 'angle_ideal'
        lines.append("    {} = {{:.{}f}}".format(keyword_ideal, digit_after_dot).format(restraint.value))
        sigma = self.dist_sigma_value(restraint.sigma) if kind == 'bond' else self.angle_sigma_value(restraint.sigma)
        lines.append("    sigma = {{:.{}f}}".format(digit_after_dot).format(sigma))
        if kind == 'bond':
            lines.append("    slack = None")
        lines.append("  }")
        return '\n'.join(lines)

    def get_dist(self, restraint, all_restraints):
        return self._prepare_restraint('bond', 3, restraint)

    def get_angle(self, restraint, all_restraints):
        return self._prepare_restraint('angle', 1, restraint)


class RefmacPrinter(PrinterBase):

    @classmethod
    def _get_alt_loc(cls, atom):
        return ' altecode {}'.format(atom.alt_loc) if atom.alt_loc != '' else ''

    @classmethod
    def _comment(cls, restraint):
        if restraint.condition_name is not None and restraint.name is not None:
            if restraint.type == 'angle':
                return '# Restraint {} {} {} {:.1f} {:.1f}'.format(
                    restraint.condition_name, restraint.type, restraint.name, restraint.value, restraint.sigma)
            return '# Restraint {} {} {} {:.3f} {:.3f}'.format(
                restraint.condition_name, restraint.type, restraint.name, restraint.value, restraint.sigma)
        return ''

    def get_dist(self, restraint, all_restraints):
        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]

        lines = []
        comment = self._comment(restraint)
        if comment != '':
            lines.append(comment)

        #print(type(restraint.value), type(restraint.sigma))
        line = "exte dist first chain {} resi {} atom {}{} second chain {} resi {} atom {}{} value {:.3f} sigma {:.3f} type {}".format(
            atom0.chain_id,
            atom0.res_id,
            atom0.atom_name,
            self._get_alt_loc(atom0),

            atom1.chain_id,
            atom1.res_id,
            atom1.atom_name,
            self._get_alt_loc(atom1),

            restraint.value,
            restraint.sigma,
            0  # type 0=override, 1=add to existing, 2=longer range distance
        )
        lines.append(line)
        return "\n".join(lines)

    def get_angle(self, restraint, all_restraints):
        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]
        atom2 = restraint.atoms[2]

        lines = []
        comment = self._comment(restraint)
        if comment != '':
            lines.append(comment)

        #print(type(restraint.value), type(restraint.sigma))
        line = "exte angle first chain {} resi {} atom {}{} next chain {} resi {} atom {}{} next chain {} resi {} atom {}{} value {:.1f} sigma {:.1f} type {}".format(
            atom0.chain_id,
            atom0.res_id,
            atom0.atom_name,
            self._get_alt_loc(atom0),

            atom1.chain_id,
            atom1.res_id,
            atom1.atom_name,
            self._get_alt_loc(atom1),

            atom2.chain_id,
            atom2.res_id,
            atom2.atom_name,
            self._get_alt_loc(atom2),

            restraint.value,
            restraint.sigma,
            0  # type 0=override, 1=add to existing, 2=longer range distance
        )
        lines.append(line)
        return "\n".join(lines)


class BusterPrinter(PrinterBase):

    @classmethod
    def _get_alt_loc(cls, atom):
        return '.{}'.format(atom.alt_loc) if atom.alt_loc != '' else ''

    @classmethod
    def _comment(cls, restraint):
        if restraint.condition_name is not None and restraint.name is not None:
            if restraint.type == 'angle':
                return '# Restraint {} {} {} {:.1f} {:.1f}'.format(
                    restraint.condition_name, restraint.type, restraint.name, restraint.value, restraint.sigma)
            return '# Restraint {} {} {} {:.3f} {:.3f}'.format(
                restraint.condition_name, restraint.type, restraint.name, restraint.value, restraint.sigma)
        return ''

    def angle_sigma_value(self, value):
        """ Default sigma for angles in BUSTER is 1.5"""
        return 1.5 if self.override_sigma else value

    def dist_sigma_value(self, value):
        """ Default sigma for bond distances in BUSTER is 0.020"""
        return 0.020 if self.override_sigma else value

    def get_dist(self, restraint, all_restraints):
        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]

        lines = []
        comment = self._comment(restraint)
        if comment != '':
            lines.append(comment)
        
        line = "NOTE BUSTER_DISTANCE ={:.3f} {:.3f} {}|{}:{}{} {}|{}:{}{}".format(
            restraint.value,
            self.dist_sigma_value(restraint.sigma),
            
            atom0.chain_id,
            atom0.res_id,
            atom0.atom_name,
            self._get_alt_loc(atom0),

            atom1.chain_id,
            atom1.res_id,
            atom1.atom_name,
            self._get_alt_loc(atom1),
        )
        lines.append(line)
        return "\n".join(lines)

    def get_angle(self, restraint, all_restraints):
        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]
        atom2 = restraint.atoms[2]

        lines = []
        comment = self._comment(restraint)
        if comment != '':
            lines.append(comment)

        line = "NOTE BUSTER_UTILANGLE {:.1f} {:.1f} {}|{}:{}{} {}|{}:{}{} {}|{}:{}{}".format(
            restraint.value,
            self.angle_sigma_value(restraint.sigma),

            atom0.chain_id,
            atom0.res_id,
            atom0.atom_name,
            self._get_alt_loc(atom0),

            atom1.chain_id,
            atom1.res_id,
            atom1.atom_name,
            self._get_alt_loc(atom1),

            atom2.chain_id,
            atom2.res_id,
            atom2.atom_name,
            self._get_alt_loc(atom2),
        )
        lines.append(line)
        return "\n".join(lines)


class CsvPrinter(PrinterBase):

    @classmethod
    def _get_alt_loc(cls, atom):
        return atom.alt_loc

    @classmethod
    def header(cls):
        return 'type,condition_name,restraint_name,chain1,resi1,atom1,altloc1,chain2,resi2,atom2,altloc2,chain3,resi3,atom3,altloc3,value,sigma'

    def get_dist(self, restraint, all_restraints):
        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]

        # "type,condition_name,restraint_name,chain1,resi1,atom1,altloc1,chain2,resi2,atom2,altloc2,chain3,resi3,atom3,altloc3,value,sigma"
        line = ",".join((str(_) for _ in (
            restraint.type,
            restraint.condition_name,
            restraint.name,

            atom0.chain_id,
            atom0.res_id,
            atom0.atom_name,
            self._get_alt_loc(atom0),

            atom1.chain_id,
            atom1.res_id,
            atom1.atom_name,
            self._get_alt_loc(atom1),

            '', '', '', '',

            '{:.3f}'.format(restraint.value),
            '{:.3f}'.format(restraint.sigma),
        )))

        return line

    def get_angle(self, restraint, all_restraints):
        atom0 = restraint.atoms[0]
        atom1 = restraint.atoms[1]
        atom2 = restraint.atoms[2]

        # "type,condition_name,restraint_name,chain1,resi1,atom1,altloc1,chain2,resi2,atom2,altloc2,chain3,resi3,atom3,altloc3,value,sigma"
        line = ",".join((str(_) for _ in (
            restraint.type,
            restraint.condition_name,
            restraint.name,

            atom0.chain_id,
            atom0.res_id,
            atom0.atom_name,
            self._get_alt_loc(atom0),

            atom1.chain_id,
            atom1.res_id,
            atom1.atom_name,
            self._get_alt_loc(atom1),

            atom2.chain_id,
            atom2.res_id,
            atom2.atom_name,
            self._get_alt_loc(atom2),

            '{:.1f}'.format(restraint.value),
            '{:.1f}'.format(restraint.sigma),
        )))

        return line
