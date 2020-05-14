#!/usr/bin/env python

import os

if config["test"]:
    config["maxSpotId"] = "-X {spots}".format(spots=config["spots"])
else:
    config["maxSpotId"] = ""