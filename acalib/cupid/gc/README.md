# PyCupid: GaussClumps

## Build
Use ```make``` to build the project with the included Makefile. This will build Starlink's ```ast``` library first, followed
by the ```pycupid``` module.

## Requirements
In order to build, the following packages need to be installed on your Python environment:

* Cython
* Numpy

## How to use
After building, you must add the path ```./libraries/ast/lib``` to the environmental variable LD_LIBRARY_PATH. 

```
$ export LD_LIBRARY_PATH=./libraries/ast/lib
```

Then, you can import the ```pycupid``` module from a Python shell 
or script. There's a ```test.py``` file included that shows how to import and calling functions from the module.
