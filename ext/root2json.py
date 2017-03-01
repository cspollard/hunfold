from sys import argv, stdout
myargv = argv[:]

import ROOT
import json

def main(args):
    dout = {}

    fjson = file(args[1])
    d = json.load(fjson)
    fjson.close()

    # get the data info
    print "loading data hist"; stdout.flush()
    datafile = ROOT.TFile.Open(d["data"]["file"])
    datahist = datafile.Get(d["data"]["hist"])
    dout["Data"] = th1ToList(datahist)

    # get the info for the nominal model
    print "loading nominal model"; stdout.flush()
    dnom = d["nominal"]
    lumi = float(dnom["lumi"])

    dout["Nominal"] = { "Lumi" : lumi, "Bkgs" : {} }

    dbkgs = dnom["backgrounds"]

    for (tbkg, dbkg) in dbkgs.iteritems():
        f = ROOT.TFile(dbkg["file"])
        h = f.Get(dbkg["hist"])
        h.Scale(1.0/lumi)
        dout["Nominal"]["Bkgs"][tbkg] = th1ToList(h)
        continue

    fmig = ROOT.TFile.Open(dnom["migration"]["file"])
    mig = th2ToList(fmig.Get(dnom["migration"]["hist"]))
    dout["Nominal"]["Mig"] = mig
    dout["Nominal"]["Sig"] = map(lambda x: 0, mig)

    # get the info for the variations
    dout["ModelVars"] = {}
    dvars = d["variations"]


    for (tvar, dvar) in dvars.iteritems():
        dtmp = {}
        if



    print json.dumps(dout); stdout.flush()

    return


def th1ToList(h):
    return [h.GetBinContent(iBin) for iBin in range(1, h.GetNbinsX()+1)]

def th2ToList(h):
    return \
        [ [ h.GetBinContent(iBin, jBin)
            for jBin in range (1, h.GetNbinsY()+1)
          ] for iBin in range (1, h.GetNbinsX()+1)
        ]

if __name__ == "__main__":
    main(myargv)
