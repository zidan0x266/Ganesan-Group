"""
Copyright (C) 2018-2020 Zidan Zhang <zhangzidan@gmail.com>
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

import numpy as np

def sigmaNerestEinstein(factor, volume, dcation, danion, dlithium, ncation, nanion, nlithium):
    """
    Nernst-Einstein conductivity
    """
    kB = 1.3806485279e-23
    T = 600
    qc = 1.602176565e-19 * factor
    qa = 1.602176565e-19 * factor
    ql = 1.602176565e-19 * factor
    number_density = (1 / (volume * 1e-30 * kB * T)) * 1e-4
    sigmaTotal = (ncation * qc**2 * dcation + nanion * qa**2 * danion + nlithium * ql**2 * dlithium) * number_density
    return sigmaTotal


def main():
    ncation = 300
    nanion = 300
    nlithium = 0
    vol = np.loadtxt('vol_c000.dat')
    dcation = np.loadtxt('Dc_c000.dat')
    danion = np.loadtxt('Da_c000.dat')
    #dlithium = np.loadtxt('Dl_c010.dat')
    dlithium = np.zeros(len(vol))
    conduct = []
    volume = [x**3 for x in vol]
    for i in range(len(vol)):
        conduct.append(sigmaNerestEinstein(0.8, volume[i], dcation[i], danion[i], dlithium[i], ncation, nanion, nlithium))
    print("{:10.6f} {:10.6f}".format(np.average(conduct), np.std(conduct)))


if __name__ == "__main__":
    main()