import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.lines as lines
import numpy as np
from sys import stdin, stdout, argv


do2D = "--2d" in argv

names = map(str.strip, stdin.readline().split(","))
xs_unsorted = np.loadtxt(stdin, delimiter=',')
print xs_unsorted

xs_sorted = np.flipud(xs_unsorted[xs_unsorted[:,0].argsort()])
print xs_sorted

xs = xs_sorted.transpose()

print xs.shape

# if we only have one param we need to add a dimension.
if len(xs.shape) == 1:
    xs.shape = (1, xs.shape[0])

npbiny = []
npbinerr = []
npbinnames = []
truthbinx = []
truthbiny = []
truthbinerr = []
normtruthbinx = []
normtruthbiny = []
normtruthbinerr = []
recobinx = []
recobiny = []
recobinerr = []

nentries = xs.shape[1]

# int() == floor() for positive numbers.
entry68 = int(nentries*0.68)
print entry68

for i in range(len(names)):
    param = xs[i]
    best68 = param[:entry68]
    name = names[i]

    hist, bins = np.histogram(param, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    best = param[0]
    med = np.median(param)

    q16 = np.min(best68)
    q84 = np.max(best68)

    print name, "mode, median, mean, low, high: %.3e, %.3e, %.3e, %.3e, %.3e" \
            % ( best, med, np.mean(param)
              , q16, q84)
    stdout.flush()

    fig, ax = plt.subplots()
    fig.suptitle(name)

    plt.bar(center, hist, align='center', width=width)
    yint = ax.get_yaxis().get_data_interval()
    plt.plot([best, best], [yint[0], yint[1]], color='red', lw=2)
    plt.plot([q16, q16], [yint[0], yint[1]], color='red', lw=2,
            linestyle="dashed")
    plt.plot([q84, q84], [yint[0], yint[1]], color='red', lw=2,
            linestyle="dashed")
    plt.savefig("%s.pdf" % name)
    plt.clf()


    plt.close()

    if name.startswith("recobin"):
        recobinx.append(float(name[7:]))
        recobiny.append(best)
        recobinerr.append((best-q16, q84-best))

    elif name.startswith("truthbin"):
        truthbinx.append(float(name[8:]))
        truthbiny.append(best)
        truthbinerr.append((best-q16, q84-best))

    elif name.startswith("normtruthbin"):
        normtruthbinx.append(float(name[12:]))
        normtruthbiny.append(best)
        normtruthbinerr.append((best-q16, q84-best))

    elif "llh" not in name:
        npbiny.append(best)
        npbinerr.append((best-q16, q84-best))
        npbinnames.append(name)

    if do2D:
        for j in range(i+1, len(names)):
            paramy = xs[j]
            namey = names[j]
            fig = plt.figure()
            fig.suptitle(name + " vs " + namey)
            plt.hist2d(param, paramy, bins=50)
            plt.colorbar()
            plt.show()
            plt.savefig("%svs%s.pdf" % (name, namey))
            plt.clf()
            plt.close()

fig = plt.figure()
fig.suptitle("reco")

plt.errorbar(recobinx, recobiny, yerr=zip(*recobinerr), xerr=0.5, fmt='o')
plt.savefig("recobin.pdf")
plt.clf()
plt.close()

fig = plt.figure()
fig.suptitle("truth")

plt.errorbar(truthbinx, truthbiny, yerr=zip(*truthbinerr), xerr=0.5, fmt='o')
plt.savefig("truthbin.pdf")
plt.clf()
plt.close()

fig = plt.figure()
fig.suptitle("normtruth")

plt.errorbar(normtruthbinx, normtruthbiny, yerr=zip(*normtruthbinerr), xerr=0.5, fmt='o')
plt.savefig("normtruthbin.pdf")
plt.clf()
plt.close()

fig = plt.figure()
fig.suptitle("nuisance params")

(npbinnames, npbiny, npbinerr) = \
        map(list, zip(*sorted(zip(npbinnames, npbiny, npbinerr))))

binsx = range(0, len(npbinnames))

ax = plt.subplots()[1]
ax.set_xticks([-1] + binsx + [len(binsx)])
ax.set_xticklabels([""] + npbinnames + [""], rotation=90,
        rotation_mode="anchor", ha="right", fontsize=8)

plt.errorbar(binsx, npbiny, yerr=zip(*npbinerr), xerr=0.5, fmt='o')
plt.tight_layout()
plt.savefig("npbin.pdf")
plt.clf()
plt.close()
