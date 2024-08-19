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
from scipy.optimize import curve_fit
from scipy.special import gamma
from scipy import stats

def relaTau(popt):
    return popt [0] * gamma(1 + 1.0 / popt[1])


def singleKWW(time, tstar, beta):
    return np.exp(-(time / tstar)**beta)


def linearKWW(time, acf):
    from scipy import stats
    X0 = [np.log(x) for x in time]
    Yt = [-np.log(x) for x in acf]
    Y0 = [np.log(x) for x in Yt]
    slope, intercept, r_value, p_value, std_err = stats.linregress(X0, Y0)
    tstar = np.exp(-intercept / slope)
    return [tstar, slope], intercept, r_value**2


def main():
    folder = ['s01', 's02', 's03', 's04', 's05',
              's06', 's07', 's08', 's09', 's10']
    tauPhi = []
    tauPhi_R2 = []
    for dir in folder:
        data = np.loadtxt('../../' + dir + '/CONCENTRATION/dihedral/acf_dih.xvg')
        params = linearKWW(data[:,0][1:] / 100, data[:,1][1:])
        tauPhi.append(relaTau(params[0]))
        tauPhi_R2.append(params[2])
    tauPhi_avg = np.average(tauPhi)
    tauPhi_std = np.std(tauPhi)
    print(tauPhi_avg, tauPhi_std)


if __name__ == "__main__":
    main()