#!/usr/bin/env python
"""
Copyright (C) 2018-2019 Jakub Krajniak <jkrajniak@gmail.com>
Zidan Zhang <zhangzidan@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import sys

__doc__ = "Tool functions."


class MyArgParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(MyArgParser, self).__init__(*args, **kwargs)

    def convert_arg_line_to_args(self, line):
        for arg in line.split():
            t = arg.strip()
            if not t:
                continue
            if t.startswith('#'):
                break
            if not t.startswith('--'):
                t = '--{}'.format(t)
            yield t

    @staticmethod
    def save_to_file(output_file, namespace):
        """Saves arguments to file so it can be read again.

        Args:
            output_file: The string with the name of output file.
            namespace: The namespace with arguements.
        """
        with open(output_file, "w") as of:
            for k, v in namespace.__dict__.iteritems():
                if v is not None:
                    of.write('{}={}\n'.format(k, v))


def _args_analysis():
    parser = MyArgParser(description='input parameters for UTAnalysis', fromfile_prefix_chars='@')
    parser.add_argument('--top', help='Input topology file')
    parser.add_argument('--trj', help='Input trajectory file')
    parser.add_argument('--begin', type=int, default=0, help='Frame to begin with (in time unit ps)')
    parser.add_argument('--end', type=int, default=0, help='Frame to end with (in time unit ps)')
    parser.add_argument('--func', type=str, default='copoly', choices=['copoly', 'homopoly', 'hopping', 'aacf', 'local', 'lasso'],
                        help='job type: copoly - block/random copolymer ' +
                             'homopoly - homo polymer ' +
                             'hopping - analysis hopping events ' +
                             'aacf - ions association auto-correlation function', required=True)
    parser.add_argument('--h5', default='pils.h5', help='HDF5 filename')
    parser.add_argument('--nt', type=int, default=1, help='Number of workers')
    parser.add_argument('--h5tag', type=str, default='htt', help='dataset name for aacf and hopping analysis')
    parser.add_argument('--trestart', type=int, default=1, help='Time shifting window for calculating ACF')
    parser.add_argument('--typeacf', type=str, default='both', help='job types: both - calculate C(t) and S(t) ' +
                                                                    'ct - calculate C(t) only ' +
                                                                    'st - calculate S(t) only')
    parser.add_argument('--cutoff', type=float, default=6.3, help='cutoff for association')
    parser.add_argument('--Icutoff', type=float, default=6.1, help='cutoff for the interface')
    parser.add_argument('--counterN', type=int, default=8, help='number of counterions per polymer chain')
    parser.add_argument('--anion', help='anion selection')
    parser.add_argument('--cation', help='cation selection')
    parser.add_argument('--beadiny', help='bead to determine the interface')
    
    return parser
    
