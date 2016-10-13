current_dir = $(shell pwd)

all: build_astlib
	python setup.py build_ext --inplace

build_astlib:
	$(shell cd ./ast && ./configure --prefix="$current_dir/libraries/ast/")
	$(MAKE) -C ./ast/
	$(MAKE) -C ./ast/ install

clean:
	rm -rf build/ libraries/
	rm -f *.so