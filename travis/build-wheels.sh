#!/bin/bash
VERSIONS=(2.7)
DIST=`uname`
if [ "$DIST" = "Linux" ];then
	PLATAFORM=manylinux1_x86_64
else
	PLATAFORM=MAC
fi

cd /workspace
for VERSION in ${VERSIONS[@]}; do
echo "Looking for py$VERSION environment ($PLATAFORM)"
 (source activate py$VERSION 2>/dev/null && echo "Environment found"  && python setup.py bdist_wheel --plat-name $PLATAFORM && source deactivate)  || (echo "Environment py$VERSION not found, creating" && conda create -y --name py$VERSION python=$VERSION numpy cython && source activate py$VERSION  &&  python setup.py bdist_wheel --plat-name $PLATAFORM && source deactivate)
  if [ "$PLATAFORM" = "MAC" ]; then
  	before="MAC.whl"
  	after="macosx_10_6_intel.macosx_10_9_intel.macosx_10_9_x86_64.macosx_10_10_intel.macosx_10_10_x86_64.whl"
  	cd setup.pydist
  	head=`ls *-TMP.whl| cut -f1-4 -d"-"`
  	mv *-TMP.whl $head-$after
  	cd ..
  fi
done
ls dist
