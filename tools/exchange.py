"""
Copyright (C) 2018-2019 Zidan Zhang <zhangzidan@gmail.com>
Jakub Krajniak <jkrajniak@gmail.com>

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

from analysis import *

def main():
    bu = rawLoad('in_htt')
    it = rawLoad('it_htt')
    bu2it = []
    it2bu = []
    for i in range(len(bu) - 1):
        bu2it.append(len(set(bu[i]) & set(it[i + 1])))
        it2bu.append(len(set(bu[i + 1]) & set(it[i])))
    print(np.average(bu2it), np.average(it2bu))


if __name__ == "__main__":
    main()