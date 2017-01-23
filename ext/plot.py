import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from sys import stdin, stdout


names = map(str.strip, stdin.readline().split(","))
xs = np.loadtxt(stdin, delimiter=',').transpose()

# if we only have one param we need to add a dimension.
if len(xs.shape) == 1:
    xs.shape = (1, xs.shape[0])

for n in range(len(names)):
    param = xs[n]
    name = names[n]
    hist, bins = np.histogram(param, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    fig = plt.figure()
    fig.suptitle(name)

    plt.bar(center, hist, align='center', width=width)
    plt.savefig("%s.png" % name)
    plt.clf()

    print name, "median, mean, low, high: %.3e, %.3e, %.3e, %.3e" \
            % ( np.median(param), np.mean(param)
              , np.percentile(param, 16), np.percentile(param, 84))
    stdout.flush()

    plt.close()
