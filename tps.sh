#!/bin/bash

PYTHONPATH=`pwd`/cmonkey:$PYTHONPATH python cmonkey/tps.py $@

