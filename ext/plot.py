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

xs_sorted = np.flipud(xs_unsorted[xs_unsorted[:,0].argsort()])

bests = xs_sorted[0]
bestfits = xs_unsorted[0]

xs = xs_sorted.transpose()
xs_unsorted = xs_unsorted.transpose()

print("parameter names:")
print(names)
print()
print("the most-likely parameter values:")
print(bests)
print()
print("the best-fit parameter values:")
print(bestfits)
print()

print("the difference between most-likely and best-fit:")
print(bests - bestfits)
print()

# if we only have one param we need to add a dimension.
if len(xs.shape) == 1:
    xs.shape = (1, xs.shape[0])

npbiny = []
npbinmed = []
npbinmederr = []
npbinerr = []
npbinnames = []
truthbinx = []
truthbiny = []
truthbinerr = []
truthbinmed = []
truthbinmederr = []
normtruthbinx = []
normtruthbiny = []
normtruthbinerr = []
normtruthbinmed = []
normtruthbinmederr = []
recobinx = []
recobiny = []
recobinerr = []
recobinmed = []
recobinmederr = []

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
    (q16, q84) = np.percentile(param, [16, 84])
    globq16 = np.min(best68)
    globq84 = np.max(best68)

    print name, "mode, median, mean, low, high: %.3e, %.3e, %.3e, %.3e, %.3e" \
            % ( best, med, np.mean(param)
              , q16, q84)
    stdout.flush()

    mean = np.mean(param)
    var = np.var(param)
    meandiff = param - mean
    autocov = np.sum(meandiff[:-1]*meandiff[1:]) / (len(meandiff)-1)
    autocorr = autocov / var

    print("stddev:")
    print(np.sqrt(var))
    print("variance:")
    print(var)
    print("autocovariance:")
    print(autocov)
    print("autocorrelation:")
    print(autocorr)
    stdout.flush()


    fig, ax = plt.subplots()
    fig.suptitle(name)

    plt.bar(center, hist, align='center', width=width, color='green')
    yint = ax.get_yaxis().get_data_interval()
    ymin = yint[0]
    ymax = yint[1]*1.2
    plt.plot([best, best], [ymin, ymax], color='blue', lw=4)
    plt.plot([med, med], [ymin, ymax], color='red', lw=2)
    plt.plot([q16, q16], [ymin, ymax], color='red', lw=2,
            linestyle="dashed")
    plt.plot([q84, q84], [ymin, ymax], color='red', lw=2,
            linestyle="dashed")
    plt.savefig("%s.pdf" % name)
    plt.clf()


    fig, axis = plt.subplots()
    fig.suptitle(name)

    plt.plot(xs_unsorted[i])
    plt.savefig("%s_t.pdf" % name)

    plt.clf()
    plt.close()


    if name.startswith("recobin"):
        recobinx.append(float(name[7:]))
        recobiny.append(best)
        recobinerr.append((best-globq16, globq84-best))
        recobinmed.append(med)
        recobinmederr.append((med-q16, q84-med))

    elif name.startswith("truthbin"):
        truthbinx.append(float(name[8:]))
        truthbiny.append(best)
        truthbinerr.append((best-globq16, globq84-best))
        truthbinmed.append(med)
        truthbinmederr.append((med-q16, q84-med))

    elif name.startswith("normtruthbin"):
        normtruthbinx.append(float(name[12:]))
        normtruthbiny.append(best)
        normtruthbinerr.append((best-globq16, globq84-best))
        normtruthbinmed.append(med)
        normtruthbinmederr.append((med-q16, q84-med))

    elif "llh" not in name:
        npbiny.append(best)
        npbinerr.append((best-globq16, globq84-best))
        npbinnames.append(name)
        npbinmed.append(med)
        npbinmederr.append((med-q16, q84-med))

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

plt.errorbar(recobinx, recobiny, xerr=0.5, fmt='o')
plt.errorbar(recobinx, recobinmed, yerr=zip(*recobinmederr), fmt='o',
        color="red")
plt.savefig("recobin.pdf")
plt.clf()
plt.close()

fig = plt.figure()
fig.suptitle("truth")

plt.errorbar(truthbinx, truthbiny, xerr=0.5, fmt='o')
plt.errorbar(truthbinx, truthbinmed, yerr=zip(*truthbinmederr),
        fmt='o', color="red")
plt.savefig("truthbin.pdf")
plt.clf()
plt.close()

fig = plt.figure()
fig.suptitle("normtruth")

plt.errorbar(normtruthbinx, normtruthbiny, xerr=0.5, fmt='o')
plt.errorbar(normtruthbinx, normtruthbinmed,
        yerr=zip(*normtruthbinmederr), fmt='o', color="red")
plt.savefig("normtruthbin.pdf")
plt.clf()
plt.close()

fig = plt.figure()
fig.suptitle("nuisance params")

(npbinnames, npbiny, npbinmed, npbinerr, npbinmederr) = \
        map(list, zip(*sorted(zip(npbinnames, npbiny, npbinmed,
            npbinerr, npbinmederr))))

nbins = len(npbinnames)
binsx = range(0, nbins)

ax = plt.subplots()[1]
ax.set_xticks([-1] + binsx + [len(binsx)])
ax.set_xticklabels([""] + npbinnames + [""], rotation=90,
        rotation_mode="anchor", ha="right", fontsize=8)

plt.errorbar(binsx, npbiny, xerr=0.5, fmt='o')
plt.errorbar(binsx, npbinmed, yerr=zip(*npbinmederr), color="red",
        fmt='o')

plt.plot([0, nbins], [-1, -1], color='gray', lw=2,
        linestyle="dashed", alpha=0.5)

plt.plot([0, nbins], [0, 0], color='gray', lw=2,
        linestyle="dashed", alpha=0.5)

plt.plot([0, nbins], [1, 1], color='gray', lw=2,
        linestyle="dashed", alpha=0.5)

plt.tight_layout()
plt.savefig("npbin.pdf")
plt.clf()
plt.close()
