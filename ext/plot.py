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

truthbinx = []
truthbiny = []
truthbinerr = []
recobinx = []
recobiny = []
recobinerr = []
for i in range(len(names)):
    param = xs[i]
    name = names[i]
    hist, bins = np.histogram(param, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    fig = plt.figure()
    fig.suptitle(name)

    plt.bar(center, hist, align='center', width=width)
    plt.savefig("%s.png" % name)
    plt.clf()

    med = np.median(param)
    q16 = np.percentile(param, 16)
    q84 = np.percentile(param, 84)

    print name, "median, mean, low, high: %.3e, %.3e, %.3e, %.3e" \
            % ( med, np.mean(param)
              , q16, q84)
    stdout.flush()

    plt.close()

    if name.startswith("recobin"):
        recobinx.append(len(recobinx)+1)
        recobiny.append(med)
        recobinerr.append((med-q16, q84-med))

    elif name.startswith("truthbin"):
        truthbinx.append(len(truthbinx)+1)
        truthbiny.append(med)
        truthbinerr.append((med-q16, q84-med))


fig = plt.figure()
fig.suptitle("reco")

plt.errorbar(recobinx, recobiny, yerr=zip(*recobinerr), xerr=0.5, fmt='o')
# plt.yscale("log")
plt.savefig("recobin.png")
plt.clf()
plt.close()

fig = plt.figure()
fig.suptitle("truth")

plt.errorbar(truthbinx, truthbiny, yerr=zip(*truthbinerr), xerr=0.5, fmt='o')
# plt.yscale("log")
plt.savefig("truthbin.png")
plt.clf()
plt.close()
