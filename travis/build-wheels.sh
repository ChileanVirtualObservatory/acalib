#!/bin/bash
set -e -x
VERSIONS=(2.7 3.4 3.5 3.6)
DIST=`uname`
if [ "$DIST" = "Linux" ];then
	PLATAFORM=manylinux1_x86_64
else
	PLATAFORM=TMP
fi

for VERSION in ${VERSIONS[@]}; do
echo "Looking for py$VERSION environment"
echo $PATH
 (source activate py$VERSION 2>/dev/null && echo "Environment found" && python /workspace/setup.py bdist_wheel --plat-name $PLATAFORM && source deactivate)  || (echo "Environment py$VERSION not found, creating" && conda create -y --name py$VERSION python=$VERSION numpy cython && source activate py$VERSION && python /workspace/setup.py bdist_wheel --plat-name $PLATAFORM && source deactivate)
  if [ "$PLATAFORM" = "TMP" ]; then
  	before="TMP.whl"
  	after="macosx_10_6_intel.macosx_10_9_intel.macosx_10_9_x86_64.macosx_10_10_intel.macosx_10_10_x86_64.whl"
  	cd /workspace/dist
  	head=`ls *-TMP.whl| cut -f1-4 -d"-"`
  	mv *-TMP.whl $head-$after
  	cd ..
  fi
done
pip install twine
