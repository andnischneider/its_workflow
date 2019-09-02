#!/usr/bin/env python

import pandas as pd


def get_samples(f, pool=False):
    samples = {}
    df = pd.read_csv(f, sep="\t", header=0)
    for i in df.index:
        r = df.loc[i]
        try:
            acc = r.accession
        except KeyError:
            acc = r.sampleName
        try:
            samples[r.sampleName][r.sampleRun] = {r.geneType: acc}
        except KeyError:
            samples[r.sampleName] = {r.sampleRun: {r.geneType: acc}}
    return samples