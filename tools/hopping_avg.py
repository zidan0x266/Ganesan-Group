import numpy as np

def outPrint(filename, xdata, ydata, yerror):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# N P(N)', file=anaout)
        for i in range(0, len(xdata)):
            print('{:10.3f} {:10.5f} {:10.5f}'.format(xdata[i], ydata[i], yerror[i]), file=anaout)


def main():
    folder = ['s01', 's02', 's03', 's04', 's05',
              's06', 's07', 's08', 's09', 's10']
    intraH, interH = [], []
    for dir in folder:
        data = np.loadtxt('../' + dir + '/CONCENTRATION/analysis_PAIR/association.dat')
        intraH.append(np.average(data[:,0]))
        interH.append(np.average(data[:,1]))
    xdata = np.linspace(0, 9, 10)
    ydata = np.average(samples, axis = 0)
    yerror = np.std(samples, axis = 0)
    outPrint('assoChain_PAIR_CONCENTRATION', xdata, ydata, yerror)


if __name__ == "__main__":
    main()