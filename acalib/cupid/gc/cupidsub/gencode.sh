## Process this file with automake to produce Makefile.in

for x in $(ls *.cgen); do 
   filename=$(echo $x|cut -d "." -f1)
	echo "/*" > $filename.c
	echo "*  Name:" >> $filename.c
	echo "*    $filename.c" >> $filename.c
	echo "" >> $filename.c
	echo "*  Purpose:" >> $filename.c
	echo "*    This file expands the generic C code held in $x to provide" >> $filename.c
	echo "*    the required type-specific implementations which can be called by" >> $filename.c
	echo "*    other functions." >> $filename.c
	echo "" >> $filename.c
	echo "*  Notes:" >> $filename.c
	echo "*    - This file is generated automatically at build time (see" >> $filename.c
	echo "*    cupidsub/Makefile.am)" >> $filename.c
	echo "*/" >> $filename.c
	echo "" >>  $filename.c
	echo "#include \"prm_par.h\"" >> $filename.c
	echo "#include \"cgeneric.h\"" >> $filename.c
	echo "" >> $filename.c
	echo "#define CGEN_CODE_TYPE CGEN_DOUBLE_TYPE" >> $filename.c
	echo "#include \"cgeneric_defs.h\"" >> $filename.c
	echo "#include \"$x\"" >> $filename.c
	echo "#undef CGEN_CODE_TYPE" >> $filename.c
	echo "" >> $filename.c
	echo "#define CGEN_CODE_TYPE CGEN_FLOAT_TYPE" >> $filename.c
	echo "#include \"cgeneric_defs.h\"" >> $filename.c
	echo "#include \"$x\"" >> $filename.c
	echo "#undef CGEN_CODE_TYPE" >> $filename.c
	echo "" >> $filename.c
done



