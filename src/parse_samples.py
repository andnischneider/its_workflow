#!/usr/bin/env python

import pandas as pd


def get_samples(f, pool=False):
    samples = {}
    df = pd.read_csv(f, sep="\t", header=0)
    genes = list(df.geneType.unique())
    for i in df.index:
        r = df.loc[i]
        if pool:
            runID = r.sampleRun
            sampleName = r.sampleName
        else:
            runID = 1
            sampleName = "{}.{}".format(r.sampleName,r.sampleRun)
        try:
            acc = r.accession
        except KeyError:
            acc = sampleName
        try:
            samples[sampleName][runID] = {r.geneType: acc}
        except KeyError:
            samples[sampleName] = {runID: {r.geneType: acc}}
    return samples, genes