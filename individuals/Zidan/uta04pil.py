import numpy as np

def outPrint(filename, xdata, ydata, yerror):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# N P(N)', file=anaout)
        for i in range(0, len(xdata)):
            print('{:10.3f} {:10.5f} {:10.5f}'.format(xdata[i], ydata[i], yerror[i]), file=anaout)


def main():
    func = 0  # 0 for chain association and 1 for pair assiciation
    folder = ['s01', 's02', 's03', 's04', 's05',
              's06', 's07', 's08', 's09', 's10']
    if func == 0:
        maxN = []
        for dir in folder:
            data = np.loadtxt('../' + dir + '/CONCENTRATION/analysis_PAIR/assoChain.dat', usecols = 1)
            maxN.append(len(data))
        samples = np.zeros((10, np.max(maxN)))
        counter = 0
        for dir in folder:
            data = np.loadtxt('../' + dir + '/CONCENTRATION/analysis_PAIR/assoChain.dat', usecols = 1)
            for i in range(len(data)):
                samples[counter][i] = data[i]
            counter += 1
        xdata = np.linspace(0, np.max(maxN) - 1, np.max(maxN))
        ydata = np.average(samples, axis = 0)
        yerror = np.std(samples, axis = 0)
        outPrint('assoChain_PAIR_CONCENTRATION', xdata, ydata, yerror)
    elif func == 1:
        maxN = []
        for dir in folder:
            data = np.loadtxt('../' + dir + '/CONCENTRATION/analysis_PAIR/assoPair.dat', usecols = 1)
            maxN.append(len(data))
        samples = np.zeros((10, np.max(maxN)))
        counter = 0
        for dir in folder:
            data = np.loadtxt('../' + dir + '/CONCENTRATION/analysis_PAIR/assoPair.dat', usecols = 1)
            for i in range(len(data)):
                samples[counter][i] = data[i]
            counter += 1
        xdata = np.linspace(0, np.max(maxN) - 1, np.max(maxN))
        ydata = np.average(samples, axis = 0)
        yerror = np.std(samples, axis = 0)
        outPrint('assoPair_PAIR_CONCENTRATION', xdata, ydata, yerror)


if __name__ == "__main__":
    main()