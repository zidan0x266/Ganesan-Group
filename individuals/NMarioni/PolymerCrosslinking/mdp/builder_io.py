"""
Copyright (C) 2014-2020 Jakub Krajniak <jkrajniak@gmail.com>
2018-2020 Zidan Zhang <zhangzidan@gmail.com>

This file is part of polymerBuilder.

lab-tools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Original code by Jakub Krajniak and Zidan Zhang
# Slight modifications so the code works as intended with the crosslinking script

import collections
import copy
import datetime
import logging
import numpy
import os
import re
import sys
import warnings

try:
    import networkx
except ImportError:
    warnings.warn("networkx not found, .get_graph() method will not be available")

__doc__ = "Set of I/O classes and functions."""

logger = logging.getLogger(__name__)

Atom = collections.namedtuple(
    'Atom', [
        'atom_id',
        'name',
        'chain_name',
        'chain_idx',
        'position'
    ])


class TopoAtom(object):
    """Atom object used in TopologyFiles."""
    atom_id = None
    atom_type = None
    chain_idx = None
    chain_name = None
    name = None
    cgnr = None
    charge = None
    mass = None
    active_site = None

    def __init__(self, atom_id=None, atom_type=None, chain_idx=None,
                 chain_name=None, name=None,
                 cgnr=None, charge=None, mass=None, active_site=None):
        self.atom_id = atom_id
        self.atom_type = atom_type
        self.chain_idx = chain_idx
        self.chain_name = chain_name
        self.name = name
        self.cgnr = cgnr
        self.charge = charge
        self.mass = mass
        self.active_site = active_site

    def __repr__(self):
        return '{} ({}): {} ({}) q={}, m={}, as={}'.format(
            self.atom_id, self.chain_idx, self.name, self.chain_name, self.charge,
            self.mass, self.active_site)


def prepare_path(file_path):
    """Prepare the file to open.

    Args:
      file_path: The file path.

    Returns:
      The path to the file.
    """

    if os.path.exists(file_path):
        file_name = os.path.basename(file_path)
        dir_name = os.path.dirname(file_path)
        if not dir_name:
            dir_name = '.'
        existing_copies = [x for x in os.listdir(dir_name) if x.startswith('_%s' % file_name)]
        if existing_copies:
            max_copy_id = max([int(x.strip('_').split('.')[-1]) for x in existing_copies])
        else:
            max_copy_id = 0
        new_file_name = '_%s.%d_' % (file_name, max_copy_id + 1)
        new_file_path = os.path.join(dir_name, new_file_name)
        print('\nFound: %s, backup on: %s\n' % (file_path, new_file_path))
        os.rename(file_path, new_file_path)

    return file_path


class TopologyFile(object):
    """Reader for GROMACS .top files.

    Args:
        file_name: The input topology file.
    """

    def __init__(self, file_name):
        self.file_name = file_name
        self.title = None
        self.atoms_updated = False
        self.new_data = {
            'bonds': {},
            'angles': {},
            'dihedrals': {},
            'improper_dihedrals': {},
            'pairs': {}
            }
        self.parsers = {}
        self.writers = {}

        # chain_name -> {chain_idx -> [at1, at2...]}
        self.chains = {}  # chain_name -> {chain_idx: [at]}
        # chain_name -> {atom_name -> [at1, at2, ...]}
        self.chain_atom_names = {}  # chain_name -> {name: [at]}
        # Chain neighbours.
        #  key: (chain_name, chain_idx),
        #  value: dict key: (chain_name, chain_idx) value: reference counter
        self.chain_neighbours = collections.defaultdict(dict)
        # Current charges of atoms
        #   key: atom id
        #   value: charge
        self.current_charges = {}

        self.atoms = {}
        self.bonds = {}
        # key: atom_id
        # value: set of atom ids that are linked to the atom
        self.bonds_def = collections.defaultdict(set)
        self.angles = {}
        self.dihedrals = {}
        self.pairs = {}
        self.improper_dihedrals = {}
        self.content = None
        self.file = None

    def init(self):
        self.__init__(self.file_name)
        logger.info('Init of topology file.')
        self.chains = {}
        self.content = None
        self.file = None
        self.atoms_updated = False
        if '__state' in self.__dict__:
            del self.__dict__['__state']
            
            
class CoordinateFile(object):
    """Coordinate file object."""
    def __init__(self, file_name):
        self.file_name = file_name
        self.title = None
        self.atoms_updated = True
        self.atoms = {}
        self.file = None
        self.data = None
        self.box = None
        self.content = None
        self.scale_factor = 1.0
        self.chains = {}
        self.fragments = collections.defaultdict(dict)
        self.id_map = {}

    def init(self):
        self.__init__(self.file_name)
        logger.info('Init of coordinate file')
        self.content = None
        self.box = None
        self.file = None
        self.atoms_updated = True
        
        
class GROFile(CoordinateFile):
    def read(self):
        """Reads the .gro file and return the atom list.

        Returns:
          The dict with atoms (key: atom_id, value: atom object).
        """

        self.file = open(self.file_name, 'r')
        if not self.content:
            self.content = self.file.readlines()

        self.title = self.content[0].replace('\r\n', '').replace('\n', '')
        number_of_atoms = int(self.content[1])

        logger.info('Reading GRO file %s', self.file_name)

        for line in self.content[2:number_of_atoms + 2]:
            chain_idx = int(line[0:5].strip())
            chain_name = line[5:10].strip()
            at_name = line[10:15].strip()
            at_id = int(line[15:20].strip())
            # Need to rescale.
            pos_x = float(line[20:28].strip()) * self.scale_factor
            pos_y = float(line[28:36].strip()) * self.scale_factor
            pos_z = float(line[36:44].strip()) * self.scale_factor

            self.atoms[at_id] = (
                Atom(
                    atom_id=at_id,
                    name=at_name,
                    chain_name=chain_name,
                    chain_idx=chain_idx,
                    position=numpy.array([pos_x, pos_y, pos_z])
                ))
            self.fragments[chain_name][at_name] = self.atoms[at_id]
            if chain_name not in self.chains:
                self.chains[chain_name] = {}
            if chain_idx not in self.chains[chain_name]:
                self.chains[chain_name][chain_idx] = {}
            self.chains[chain_name][chain_idx][at_name] = self.atoms[at_id]

        # Reads the box size, the last line.
        self.box = numpy.array(
            list(map(float, [_f for _f in self.content[number_of_atoms + 2].split(' ') if _f]))
            ) * self.scale_factor

    def remove_atom(self, atom_id, renumber=True):
        """Remove atom and renumber the file."""
        atom_to_remove = self.atoms[atom_id]
        try:
            del self.fragments[atom_to_remove.chain_name][atom_to_remove.name]
            del self.chains[atom_to_remove.chain_name][atom_to_remove.chain_idx][atom_to_remove.name]
        except KeyError:
            pass
        del self.atoms[atom_id]

        if renumber:
            new_at_id = 1
            new_atoms = {}
            for at_id in self.atoms:
                new_atoms[new_at_id] = self.atoms[at_id]._replace(atom_id=new_at_id)
                new_at_id += 1
            self.atoms = new_atoms

    def renumber(self):
        """Renumber atoms with new id"""
        new_at_id = 1
        new_atoms = {}
        for at_id in self.atoms:
            new_atoms[new_at_id] = self.atoms[at_id]._replace(atom_id=new_at_id)
            new_at_id += 1
        self.atoms = new_atoms


    @staticmethod
    def copy(input_gro, particle_ids=None, renumber=False):
        """Make copy of GROFile."""
        output_gro = GROFile(input_gro.file_name)
        output_gro.box = input_gro.box
        output_gro.title = input_gro.title
        output_gro.id_map = {}
        if particle_ids:
            if renumber:
                new_pid = 1
                for pid in particle_ids:
                    at = input_gro.atoms[pid]
                    output_gro.id_map[new_pid] = pid
                    output_gro.atoms[new_pid] = Atom(
                        atom_id=new_pid,
                        name=at.name,
                        chain_name=at.chain_name,
                        chain_idx=at.chain_idx,
                        position=at.position
                    )
                    new_pid += 1
            else:
                for pid in particle_ids:
                    output_gro.atoms[pid] = copy.copy(input_gro.atoms[pid])
                    output_gro.id_map[pid] = pid
        else:
            output_gro.atoms = copy.copy(input_gro.atoms)
        return output_gro

    def write(self, file_name=None, force=False, append=False):
        """Writes the content to the output file.

        Args:
          file_name: The new file name, otherwise the old one will be used.
          force: Force to save even if any atoms were not updated.
          append: If set to true then the frame will be added to previouse one in the file.
        """

        if self.atoms_updated or force:
            output = []
            if self.title:
                output.append(self.title)
            else:
                output.append('XXX of molecules')
            # Puts the number of atoms
            output.append('%d' % len(self.atoms))
            # Puts the definition of the atoms, fixed format.
            fmt = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
            for at_id in sorted(self.atoms):
                at = self.atoms[at_id]
                output.append(fmt % (
                    at.chain_idx,
                    at.chain_name,
                    at.name,
                    at.atom_id,
                    at.position[0],
                    at.position[1],
                    at.position[2]
                    ))

            output.append('%f %f %f\n' % tuple(self.box))
            if append:
                write_file_path = file_name if file_name else self.file_name
            else:
                write_file_path = prepare_path(file_name if file_name else self.file_name)
            logger.info('Writing GRO file %s', write_file_path)
            output_file = open(write_file_path, 'a+' if append else 'w')
            output_file.writelines('\n'.join(output))
            if not append:
                output_file.write('\n')
            output_file.close()
            self.atoms_updated = False

    def update_positions(self, system, use_id_map=False):
        """Update positions."""
        for at_id in self.atoms:
            p = system.storage.getParticle(self.id_map[at_id])
            old_atom = self.atoms[at_id]
            self.atoms[at_id] = old_atom._replace(position=p.pos)

    def dump(self, system, filename, particle_ids, chain_name, chain_idx, atom_name):
        """Dump data from storage."""
        for at_id in particle_ids:
            p = system.storage.getParticle(at_id)
            self.atoms[at_id] = Atom(
                atom_id=at_id,
                name=atom_name[at_id],
                chain_name=chain_name[at_id],
                chain_idx=chain_idx[at_id],
                position=p.pos
            )
        self.title = 'XXX'
        self.box = system.bc.boxL
        self.write(filename, force=True)


class GroTopFile(object):
    """Very basic representation of topology file."""
    def __init__(self, file_name):
        self.file_name = file_name
        self.title = None
        self.atoms_updated = False
        self.new_data = {
            'bonds': {},
            'angles': {},
            'dihedrals': {},
            'improper_dihedrals': {},
            'pairs': {}
        }
        self.chains = {}
        self.chains_atoms = {}
        self.chain_atom_names = {}
        self.chain_neighbours = collections.defaultdict(dict)
        self.current_charges = {}

        self.atoms = {}
        self.bonds = {}
        # key: atom_id
        # value: set of atom ids that are linked to the atom
        self.bonds_def = collections.defaultdict(set)
        self.angles = {}
        self.dihedrals = {}
        self.pairs = {}
        self.improper_dihedrals = {}
        self.content = None
        self.file = None
        self.parsers = {
            'atoms': self._parse_atoms,
            'bonds': self._parse_bonds,
            'dihedrals': self._parse_dihedrals,
            'improper_dihedrals': self._parse_improper_dihedrals,
            'pairs': self._parse_pairs,
            'pairtypes': self._parse_pairtypes,
            'angles': self._parse_angles,
            'moleculetype': self._parse_moleculetype,
            'system': self._parse_system,
            'molecules': self._parse_molecules,
            'atomtypes': self._parse_atomtypes,
            'nonbond_params': self._parse_nonbond_params,
            'bondtypes': self._parse_bondtypes,
            'angletypes': self._parse_angletypes,
            'dihedraltypes': self._parse_dihedraltypes,
            'defaults': self._parse_defaults
        }

        self.writers = {
            'defaults': self._write_defaults,
            'moleculetype': self._write_moleculetype,
            'system': self._write_system,
            'molecules': self._write_molecules,
            'atomtypes': self._write_atomtypes,
            'bondtypes': self._write_bondtypes,
            'angletypes': self._write_angletypes,
            'dihedraltypes': self._write_dihedraltypes,
            'atoms': self._write_atoms,
            'bonds': self._write_bonds,
            'angles': self._write_angles,
            'dihedrals': self._write_dihedrals,
            'improper_dihedrals': self._write_improper_dihedrals,
            'pairs': self._write_pairs,
            }
        self.current_charges = {}
        self.atomtypes = {}
        self.nonbond_params = {}
        self.bondtypes = {}
        self.angletypes = {}
        self.dihedraltypes = {}
        self.pairtypes = {}
        self.header_section = []
        self.defaults = None
        self.moleculetype = []
        self.molecules = []
        self.system_name = None

    def init(self):
        """Reset the class properties without creating the object again."""
        logger.info('Init of topology file.')
        self.new_data = {
            'bonds': {},
            'angles': {},
            'dihedrals': {},
            'improper_dihedrals': {},
            'pairs': {}
        }

        if '__state' in self.__dict__:
            del self.__dict__['__state']
        self.current_charges = {}
        self.atoms_updated = False

    def get_graph(self):
        """Returns graph."""
        output_graph = networkx.Graph(box=None)
        for at_id, g_at in self.atoms.items():
            output_graph.add_node(
                at_id,
                name=g_at.name,
                res_id=g_at.chain_idx,
                position=(-1, -1, -1),
                chain_name=g_at.chain_name)

        for (b1, b2), params in self.bonds.items():
            output_graph.add_edge(b1, b2, params=params)

        if 'bonds' in self.new_data:
            for (b1, b2), params in self.new_data['bonds'].items():
                output_graph.add_edge(b1, b2, params=params)

        for n_id in output_graph.nodes():
            output_graph.node[n_id]['degree'] = output_graph.degree(n_id)
        return output_graph

    def update_position(self, pdbfile):
        """Reads the position data from the coordinate file and update the atoms.

        Args:
          pdbfile: The pdb file.
        """
        logger.info('Update position from file %s', pdbfile.file_name)
        for k, v in pdbfile.atoms.items():
            self.atoms[k].position = v.position

    def remove_atom(self, atom_id, renumber=True):
        """Removes atom from topology and clean data structures."""

        atom_to_remove = self.atoms[atom_id]
        try:
            self.chains[atom_to_remove.chain_name][atom_to_remove.chain_idx].remove(atom_to_remove)
            self.chain_atom_names[atom_to_remove.chain_name][atom_to_remove.name].remove(atom_to_remove)
        except KeyError:
            pass
        del self.atoms[atom_id]

        # Renumber data.
        old2new = {k: k for k in self.atoms}
        if renumber:
            new_at_id = 1
            old2new = {}
            new_atoms = {}
            for at_id in sorted(self.atoms):
                new_atoms[new_at_id] = self.atoms[at_id]
                new_atoms[new_at_id].atom_id = new_at_id
                old2new[at_id] = new_at_id
                new_at_id += 1
            self.atoms = new_atoms

        # Clean bonded structures.
        self.bonds = {k: v for k, v in list(self.bonds.items()) if atom_id not in k}
        self.angles = {k: v for k, v in list(self.angles.items()) if atom_id not in k}
        self.dihedrals = {k: v for k, v in list(self.dihedrals.items()) if atom_id not in k}
        self.pairs = {k: v for k, v in list(self.pairs.items()) if atom_id not in k}
        self.improper_dihedrals = {k: v for k, v in list(self.improper_dihedrals.items()) if atom_id not in k}

        # And new_data
        for k in self.new_data:
            self.new_data[k] = {p: v for p, v in list(self.new_data[k].items()) if atom_id not in p}


    def renumber(self):
        """Renumber topology"""
        # Clean bonded structures.
        old2new = {}
        new_at_id = 1
        new_atoms = {}
        for at_id in sorted(self.atoms):
            new_atoms[new_at_id] = self.atoms[at_id]
            new_atoms[new_at_id].atom_id = new_at_id
            old2new[at_id] = new_at_id
            new_at_id += 1
        self.atoms = new_atoms

        #self.bonds = {tuple(map(old2new.get, k)): v
        #              for k, v in list(self.bonds.items())}
        #self.angles = {tuple(map(old2new.get, k)): v
        #               for k, v in list(self.angles.items())}
        #self.dihedrals = {tuple(map(old2new.get, k)): v
        #                  for k, v in list(self.dihedrals.items())}
        #self.pairs = {tuple(map(old2new.get, k)): v
        #              for k, v in list(self.pairs.items())}
        #self.improper_dihedrals = {tuple(map(old2new.get, k)): v
        #                           for k, v in list(self.improper_dihedrals.items())}
        #
        ## And new_data
        #for k in self.new_data:
        #    self.new_data[k] = {tuple(map(old2new.get, p)): v for p, v in list(self.new_data[k].items())}

    def _replicate_lists(self, n_mols, n_atoms, input_list, shift=0):
        return {
            tuple([shift+x+(mol*n_atoms) for x in l]): v
            for mol in range(n_mols) for l, v in list(input_list.items())
            }

    def replicate(self):
        """Replicate molecules"""
        nmols = int(self.molecules[0]['mol'])
        atoms = copy.copy(self.atoms)
        for mol_id in range(1, nmols):
            for at_id in atoms:
                new_at_id = mol_id*len(atoms) + at_id
                self.atoms[new_at_id] = copy.copy(atoms[at_id])
                self.atoms[new_at_id].atom_id = new_at_id
                self.atoms[new_at_id].chain_idx = mol_id + 1

        self.bonds = self._replicate_lists(nmols, len(atoms), self.bonds.copy())
        self.angles = self._replicate_lists(nmols, len(atoms), self.angles.copy())
        self.dihedrals = self._replicate_lists(nmols, len(atoms), self.dihedrals.copy())
        self.improper_dihedrals = self._replicate_lists(nmols, len(atoms), self.improper_dihedrals.copy())

    def read(self):
        """Reads the topology file."""

        if not self.content:
            self.file = open(self.file_name, 'r')
            self.content = self.file.readlines()

        logger.info('Reading top file %s', self.file_name)

        # New version
        current_parser = None
        visited_sections = set()
        section_name = None
        previous_section = None
        new_content = []
        for line in self.content:
            line = line.strip()
            if 'include' in line:
                self.header_section.append('{}\n'.format(line.strip()))
            elif line.startswith(';') or line.startswith('#') or len(line) == 0:
                continue
            elif line.startswith('['):  # Section
                previous_section = section_name
                section_name = line.replace('[', '').replace(']', '').strip()
                # Hack for GROMACS improper_dihedrals
                if previous_section == 'dihedrals' and section_name == 'dihedrals':
                    section_name = 'improper_dihedrals'
                current_parser = self.parsers.get(section_name)
                if current_parser is not None:
                    print(('{}: Reading section {}'.format(self.file_name, section_name)))
                else:
                    print(('Parser for section {} not defined'.format(section_name)))
                visited_sections.add(previous_section)
            else:
                if current_parser is not None and section_name not in visited_sections:
                    raw_data = [_f for _f in line.split() if _f]
                    if raw_data:
                        current_parser(raw_data)  # pylint:disable=E1102

    def write(self, filename=None, force=False):
        """Updates the topology file.

        Args:
          filename: The optional output filename.
        """
        if filename is None:
            filename = self.file_name
        output_file = open(prepare_path(filename), 'w')

        new_data = []
        current_section = None
        previous_section = None
        skip_lines = False
        section_writer = None

        if self.content is None or force:
            sections = []
            if self.defaults:
                sections.append('defaults')
            if self.atomtypes or self.new_data.get('atomtypes'):
                sections.append('atomtypes')
            if self.bondtypes or self.new_data.get('bondtypes'):
                sections.append('bondtypes')
            if self.angletypes or self.new_data.get('angletypes'):
                sections.append('angletypes')
            if self.dihedraltypes or self.new_data.get('dihedraltypes'):
                sections.append('dihedraltypes')
            sections.extend([
                'moleculetype',
                'atoms',
                'bonds',
                'angles',
                'dihedrals',
                'dihedrals',
                'pairs'])
            sections.extend([
                'system',
                'molecules'
            ])
            self.content = []
            for s in sections:
                self.content.append('[ %s ]\n' % s)
                self.content.append('\n')

        new_data.extend(self.header_section)

        for line in self.content:
            tmp_line = line.strip()
            if tmp_line.startswith('['):  # section part
                new_data.append(line)
                previous_section = current_section
                current_section = tmp_line.replace('[', '').replace(']', '').strip()
                if previous_section == 'dihedrals' and current_section == 'dihedrals':
                    current_section = 'improper_dihedrals'
                section_writer = self.writers.get(current_section)
                print(('{}: Writing section {}'.format(filename, current_section)))
                skip_lines = False
            elif tmp_line.startswith(';'):
                new_data.append(line)
            elif tmp_line.startswith('#'):
                continue
            else:
                if section_writer is None:  # there is no special writer, simply copy the line
                    new_data.append(line)
                elif not skip_lines:
                    output_writer = section_writer()
                    if output_writer:
                        new_data.extend(['%s\n' % x for x in output_writer])
                    new_data.extend(['\n'])
                    skip_lines = True

        logger.info('Writing topology file %s...', filename)
        output_file.writelines(new_data)
        output_file.close()
        self.atoms_updated = False

    # Parsers for the data.
    def _parse_bonds(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:2]))
        self.bonds[atom_tuple] = raw_data[2:]

        self.bonds_def[atom_tuple[0]].add(atom_tuple[1])
        self.bonds_def[atom_tuple[1]].add(atom_tuple[0])

    def _parse_atomtypes(self, raw_data):
        #  EPO_C1 C 6 14.027000 A 0.3850 0.58576
        raw_line = ' '.join(raw_data)
        raw_line = raw_line.split(';')[0]
        before_ptype, at_type, after_ptype = re.split('\s+(A|S|D|V)+\s+', raw_line)

        before_ptype = before_ptype.split()
        name = before_ptype.pop(0)
        bonded_type = None
        atomic_number = None
        if len(before_ptype) == 2:
            mass = before_ptype.pop(0)
            charge = before_ptype.pop(0)
        else:
            if before_ptype[0].isalpha():
                bonded_type = before_ptype.pop(0)
            elif before_ptype[0].isdigit():
                atomic_number = before_ptype.pop(0)
            if len(before_ptype) == 2:
                mass = before_ptype.pop(0)
                charge = before_ptype.pop(0)
            elif len(before_ptype) == 3:
                atomic_number = before_ptype.pop(0)
                mass = before_ptype.pop(0)
                charge = before_ptype.pop(0)

        ts = after_ptype.split()
        sigma, epsilon = ts

        if at_type not in ['A', 'S', 'D', 'V']:
            print(('Wrong particle type {} in atomtypes line: {}'.format(at_type, raw_data)))

        self.atomtypes[name] = {
            'name': name,
            'mass': float(mass),
            'charge': float(charge),
            'type': at_type,
            'bonded_type': bonded_type,
            'atomic_number': atomic_number,
            'sigma': float(sigma),
            'epsilon': float(epsilon)
        }

    def _parse_pairtypes(self, raw_data):
        i, j = raw_data[:2]
        if i not in self.pairtypes:
            self.pairtypes[i] = {}
        if j not in self.pairtypes:
            self.pairtypes[j] = {}

        self.pairtypes[i][j] = {
            'func': int(raw_data[2]),
            'params': raw_data[3:]
        }

    def _parse_nonbond_params(self, raw_data):
        i, j = raw_data[:2]
        k = tuple(sorted(raw_data[:2]))
        if k in self.nonbond_params:
            raise RuntimeError('{} already exists, wrong topology'.format(k))
        self.nonbond_params[k] = {
            'func': int(raw_data[2]),
            'params': raw_data[3:]}

    def _parse_bondtypes(self, raw_data):
        i, j = raw_data[:2]
        if i not in self.bondtypes:
            self.bondtypes[i] = {}
        if j not in self.bondtypes:
            self.bondtypes[j] = {}

        params = {
            'func': int(raw_data[2]),
            'params': raw_data[3:]
        }
        self.bondtypes[i][j] = params
        self.bondtypes[j][i] = params
        self.bondtypes[(i, j)] = self.bondtypes[i][j]
        self.bondtypes[(j, i)] = self.bondtypes[i][j]

    def _parse_angletypes(self, raw_data):
        i, j, k = raw_data[:3]
        if i not in self.angletypes:
            self.angletypes[i] = {}
        if j not in self.angletypes[i]:
            self.angletypes[i][j] = {}
        if k not in self.angletypes:
            self.angletypes[k] = {}
        if j not in self.angletypes[k]:
            self.angletypes[k][j] = {}

        params = {
            'func': int(raw_data[3]),
            'params': raw_data[4:]
        }
        self.angletypes[i][j][k] = params
        self.angletypes[k][j][i] = params

        self.angletypes[(i, j, k)] = params
        self.angletypes[(k, j, i)] = params

    def _parse_dihedraltypes(self, raw_data):
        i, j, k, l = raw_data[:4]
        if i not in self.dihedraltypes:
            self.dihedraltypes[i] = {}
        if j not in self.dihedraltypes[i]:
            self.dihedraltypes[i][j] = {}
        if k not in self.dihedraltypes[i][j]:
            self.dihedraltypes[i][j][k] = {}
        if l not in self.dihedraltypes:
            self.dihedraltypes[l] = {}
        if k not in self.dihedraltypes[l]:
            self.dihedraltypes[l][k] = {}
        if j not in self.dihedraltypes[l][k]:
            self.dihedraltypes[l][k][j] = {}

        params = {
            'func': int(raw_data[4]),
            'params': raw_data[5:]
        }
        self.dihedraltypes[i][j][k][l] = params
        self.dihedraltypes[l][k][j][i] = params
        self.dihedraltypes[(l, k, j, i)] = params
        self.dihedraltypes[(i, j, k, l)] = params

    def _parse_atoms(self, raw_data):
        at = TopoAtom()
        at.atom_id = int(raw_data[0])
        at.atom_type = raw_data[1]
        at.chain_idx = int(raw_data[2])
        at.chain_name = raw_data[3]
        at.name = raw_data[4]
        at.cgnr = int(raw_data[5])

        if len(raw_data) > 6:
            at.charge = float(raw_data[6])
        if len(raw_data) > 7:
            at.mass = float(raw_data[7])

        if at.chain_name not in self.chains:
            self.chains[at.chain_name] = collections.defaultdict(list)
            self.chain_atom_names[at.chain_name] = collections.defaultdict(list)

        self.chains[at.chain_name][at.chain_idx].append(at)
        self.chain_atom_names[at.chain_name][at.name].append(at)
        self.atoms[at.atom_id] = at

    def _parse_angles(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:3]))
        self.angles[atom_tuple] = raw_data[3:]

    def _parse_dihedrals(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:4]))
        self.dihedrals[atom_tuple] = raw_data[4:]

    def _parse_improper_dihedrals(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:4]))
        self.improper_dihedrals[atom_tuple] = raw_data[4:]

    def _parse_pairs(self, raw_data):
        atom_tuple = tuple(map(int, raw_data[0:2]))
        self.pairs[atom_tuple] = raw_data[2:]

    def _parse_moleculetype(self, raw_data):
        self.moleculetype.append({
            'name': raw_data[0],
            'nrexcl': raw_data[1]
        })

    def _parse_molecules(self, raw_data):
        self.molecules.append({
            'name': raw_data[0],
            'mol': raw_data[1]
        })

    def _parse_system(self, raw_data):
        self.system_name = raw_data[0]

    def _parse_defaults(self, raw_data):
        self.defaults = {
            'nbfunc': int(raw_data[0]),
            'combinationrule': int(raw_data[1]),
            'comb-rule': int(raw_data[1])
            }
        if len(raw_data) > 2:
            self.defaults['gen-pairs'] = raw_data[2] == 'yes'
            self.defaults['fudgeLJ'] = float(raw_data[3])
            self.defaults['fudgeQQ'] = float(raw_data[4])
        else:
            self.defaults['gen-pairs'] = False
            self.defaults['fudgeLJ'] = 1.0
            self.defaults['fudgeQQ'] = 1.0

    # Writers
    def _write_atoms(self):
        return_data = []
        for atom_id in sorted(self.atoms):
            x = self.atoms[atom_id]
            return_data.append('%s %s %s %s %s %s %s %s' % (
                x.atom_id,
                x.atom_type,
                x.chain_idx,
                x.chain_name,
                x.name,
                x.cgnr,
                x.charge if x.charge is not None else '0.0',
                x.mass if x.mass is not None else ''
            ))
        return return_data

    def _write_atomtypes(self):
        return_data = []
        for atom_type, values in self.atomtypes.items():
            return_data.append('{name} {mass} {charge} {type} {sigma} {epsilon}'.format(
                **values))
        return return_data

    def _write_bondtypes(self):
        return_data = []
        for i in self.bondtypes:
            if not isinstance(i, tuple):
                for j, params in list(self.bondtypes[i].items()):
                    return_data.append('{} {} {} {}'.format(i, j, params['func'], ' '.join(params['params'])))
        return return_data

    def _write_angletypes(self):
        return_data = []
        for i in self.angletypes:
            if not isinstance(i, tuple):
                for j in self.angletypes[i]:
                    for k, params in list(self.angletypes[i][j].items()):
                        return_data.append('{} {} {} {} {}'.format(
                            i, j, k, params['func'], ' '.join(params['params'])))
        return return_data

    def _write_dihedraltypes(self):
        return_data = []
        for i in self.dihedraltypes:
            if not isinstance(i, tuple):
                for j in self.dihedraltypes[i]:
                    for k in self.dihedraltypes[i][j]:
                        for l, params in list(self.dihedraltypes[i][j][k].items()):
                            return_data.append('{} {} {} {} {} {}'.format(
                                i, j, k, l, params['func'], ' '.join(params['params'])))
        return return_data

    def _write_bonds(self):  # pylint:disable=R0201
        return_data = []
        return_data.extend(self._write_default(self.bonds))
        return_data.extend(self._write_default(self.new_data['bonds'], self.bonds))
        return return_data

    def _write_pairs(self):  # pylint:disable=R0201
        return_data = []
        return_data.extend(self._write_default(self.pairs))
        return_data.extend(
            self._write_default(self.new_data['pairs'], self.pairs)
        )
        return return_data

    def _write_angles(self):
        return_data = []
        return_data.extend(self._write_default(self.angles))
        return_data.extend(
            self._write_default(self.new_data['angles'], self.angles)
        )
        return return_data

    def _write_dihedrals(self):
        return_data = []
        return_data.extend(self._write_default(self.dihedrals))
        return_data.extend(
            self._write_default(self.new_data['dihedrals'], self.dihedrals)
        )
        return return_data

    def _write_improper_dihedrals(self):
        return_data = []
        return_data.extend(self._write_default(self.improper_dihedrals))
        return_data.extend(
            self._write_default(self.new_data['improper_dihedrals'], self.improper_dihedrals)
        )
        return return_data

    def _write_defaults(self):
        if self.defaults:
            if 'gen-pairs' in self.defaults:
                self.defaults['gen-pairs'] = 'yes' if self.defaults['gen-pairs'] else 'no'
            return [
                '{nbfunc} {comb-rule} {gen-pairs} {fudgeLJ} {fudgeQQ}'.format(**self.defaults)
            ]
        return []

    def _write_moleculetype(self):
        return ['{name} {nrexcl}\n'.format(**x) for x in self.moleculetype]

    def _write_system(self):
        return [self.system_name]

    def _write_molecules(self):
        return ['{name} {mol}\n'.format(**x) for x in self.molecules]

    def _write_default(self, datas=None, check_in=None):  # pylint:disable=R0201
        if check_in is None:
            check_in = []

        if datas is None:
            return False

        if not isinstance(datas, list):
            datas = [datas]

        if None in datas:
            return False

        flat_data = []
        for data in datas:
            for key, values in data.items():
                rev_key = tuple(reversed(key))
                if tuple(key) not in check_in or rev_key not in check_in or rev_key not in data:
                    flat_data.append(list(key) + list(values))
        flat_data.sort()
        return ['%s' % ' '.join(map(str, x)) for x in flat_data]


def read_coordinates(file_name):
    file_suffix_class = {
        'pdb': PDBFile,
        'gro': GROFile
    }
    file_suffix = file_name.split('.')[-1]
    f = file_suffix_class[file_suffix](file_name)
    f.read()
    return f


def read_topology(file_name, settings):
    file_suffix_class = {
        'top': GROMACSTopologyFile,
    }
    file_suffix = file_name.split('.')[-1]
    return file_suffix_class[file_suffix](file_name, settings)