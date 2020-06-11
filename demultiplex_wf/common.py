#!/usr/bin/env python

import os
import urllib
import json


def generate_fastq_ftp(sra_id):
    url_base = "https://www.ebi.ac.uk/ena/portal/api/filereport?"
    query = "accession={sra_id}&result=read_run".format(sra_id=sra_id)
    fmt = "&fields=fastq_ftp&format=json&download=true"
    url = "{}{}{}".format(url_base, query, fmt)
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())[0]
    return data["fastq_ftp"].split(";")


if config["test"]:
    config["maxSpotId"] = "-X {spots}".format(spots=config["spots"])
    ruleorder: fastq_dump > download_data
else:
    config["maxSpotId"] = ""
    ruleorder: download_data > fastq_dump
    config["ftp"] = {}
    sys.stderr.write("Generating ftp urls...\n")
    for sample in config["samples"].keys():
        config["ftp"][sample] = {"seqs": {"R1": "", "R2": ""},
                                 "index": {"R1": "", "R2": ""}}
        for key in ["seqs", "index"]:
            sys.stderr.write("... {} ({})\n".format(sample, key))
            r1, r2 = generate_fastq_ftp(config["samples"][sample][key])
            config["ftp"][sample][key]["1"] = r1
            config["ftp"][sample][key]["2"] = r2
