#!/bin/bash
set -e -x
VERSIONS=(cp27-cp27mu cp34-cp34m cp35-cp35m cp36-cp36m)

yum update -y
yum groupinstall -y "Development Tools"
yum install -y csh libXext-devel libXau-devel libX11-devel libXt-devel libxml2-devel ncurses-devel texlive-multirow python-devel Cython
 
# Compile wheels
for PYBIN in ${VERSIONS[@]}; do
    PYBIN=/opt/python/$PYBIN/bin/
    "${PYBIN}/pip" install -r /io/requirements.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

for PYBIN in ${VERSIONS[@]}; do
    PYBIN=/opt/python/$PYBIN/bin/
    "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/python" /io/testing/run_test.py)
done
