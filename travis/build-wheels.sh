#!/bin/bash
set -e -x

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r /workspace/requirements.txt
    "${PYBIN}/pip" wheel /workspace/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /workspace/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install python-manylinux-demo --no-index -f /workspace/wheelhouse
    (cd "$HOME"; "${PYBIN}/python" /workspace/testing/run_test.py)
done
