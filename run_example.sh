#!/usr/bin/env bash
set -e -x

snakemake exploration_all integration_all \
  --profile .profiles/local \
  --configfile \
    configs/integration/config.yaml \
    $@
