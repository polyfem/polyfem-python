#!/usr/bin/env bash
set -eux

export CMAKE_GENERATOR=Ninja
export CONDA_PREFIX="${PREFIX}"
"${PYTHON}" -m pip install . --no-deps --no-build-isolation -v
