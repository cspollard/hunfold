from array import array
from sys import stdin, stdout, argv
args = argv[:]
import ROOT

ROOT.TH1.SetDefaultSumw2(True)

names = map(str.strip, stdin.readline().split(","))

f = ROOT.TFile.Open(args[1], "CREATE")

if not f:
    print "unable to open file", f.GetName()
    print "perhaps it already exists?"
    exit(-1)

t = ROOT.TTree("stattree", "stattree")

d = []
for name in names:
    x = array('d', [0])
    d.append(x)
    t.Branch(name, x, name + "/D")


for line in stdin:
    vals = map(str.strip, line.split(','))

    for (val, x) in zip(vals, d):
        x[0] = float(val)

    t.Fill()

t.Write()

f.Close()

# hists = []
# for name in names:
#     t.Draw(name + ">>h" + name)
#     hist = f.Get("h" + name)
#     hists.append(hist)
# 
# map(ROOT.TObject.Write, hists + [t])
