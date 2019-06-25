INSTALL_DIR=/usr/local/lib/
RM= /bin/rm -vf
CC=gcc
PYTHON=python
ARCH=UNDEFINED
PWD=pwd
CDR=$(shell pwd)

EDCFLAGS:= $(CFLAGS)
EDLDFLAGS:= $(LDFLAGS)

ifeq ($(ARCH),UNDEFINED)
	ARCH=$(shell uname -m)
endif

OS=$(shell uname -s)
ifeq ($(OS),Darwin)
	LIBEXT= dylib
	EDCFLAGS:= -arch $(ARCH) $(EDCFLAGS)
	EDLDFLAGS:= -arch $(ARCH) $(EDLDFLAGS)
	LINKOPTIONS:= -dynamiclib -single_module
	RMDIR=rmdir
	ECHO=echo
	TRUE=TRUE
else
	LIBEXT= so
	LINKOPTIONS:= -shared
	RMDIR=rmdir -v
	ECHO=/bin/echo -e
	TRUE=true
endif
TARGETLIB= libmagfieldgen.$(LIBEXT)
LINKOPTIONS:=-L lib/

lib_magfieldgen_objects = src/magfieldgen.o
magfieldgen_example = src/magfieldgen.c

all: dcmt build/$(TARGETLIB) example/example

dcmt:
	cd lib && make && cd ..

build:
	mkdir build

build/$(TARGETLIB): $(lib_magfieldgen_objects) build
	$(CC) $(LINKOPTIONS) -o $@ \
	 $(EDLDFLAGS)\
	 $(lib_magfieldgen_objects) -lm -lfftw3 -fopenmp -ldcmt

example/example: $(magfieldgen_example)
	 	$(CC) -o $@ -I include/ -L lib/ \
	 	 $(EDLDFLAGS)\
	 	 $(magfieldgen_example) -lm -lfftw3 -fopenmp -ldcmt

%.o: %.c
	$(CC) -O3 -fPIC -Wall -fopenmp -c $< -o $@ -I src/ -I include/


install: build/$(TARGETLIB)
	cp $< $(INSTALL_DIR)$(TARGETLIB)

#example:
#	$(CC) $(LINKOPTIONS) -o example/example \
#	magfieldgen.c -lm -lfftw3 -fopenmp -ldcmt
#	$(shell $CDR/example/example)
# INSTALL THE PYTHON WRAPPER
pywrapper:
	(ls $(INSTALL_DIR)$(TARGETLIB) || ($(ECHO) "Cannot find library '$(INSTALL_DIR)$(TARGETLIB)', run again with 'INSTALL_DIR=' set to the directory you installed the library in" && exit -1))
	(ls py/magfieldgen/.magfieldgen.py || mv py/magfieldgen/magfieldgen.py py/magfieldgen/.magfieldgen.py)
	sed "s#TEMPLATE_LIBRARY_PATH#'$(INSTALL_DIR)'#g" py/magfieldgen_TEMPLATE.py > py/magfieldgen/magfieldgen.py
	((cd py/magfieldgen && $(ECHO) 'import magfieldgen' | $(PYTHON)) && $(ECHO) 'Successfully installed Python wrapper' || ($(ECHO) 'Something went wrong installing Python wrapper' && exit -1))
	(cd py && $(PYTHON) setup.py install --record ../python_install.txt)


#
# TEST THE INSTALLATION
#

testpy:
	(cd py/magfieldgen && $(PYTHON) magfieldgen.py)

testc:
	(cd example && time ./example)


.PHONY: clean spotless rmbuild

clean:
	$(RM) $(lib_magfieldgen_objects)
	$(RM) example/example
	cd lib && make clean && cd ..

spotless: clean rmbuild pyclean
	$(RM) src/*.o
	$(RM) py/*.pyc
	$(RM) py/magfieldgen/*.pyc
	$(RM) example/example
	(ls py/magfieldgen/.magfieldgen.py && mv py/magfieldgen/.magfieldgen.py py/magfieldgen/magfieldgen.py || $(TRUE))

pyclean:
	$(RM) -R py/build
	$(RM) -R py/dist
	$(RM) -R py/magfieldgen.egg-info
	(ls py/magfieldgen/.magfieldgen.py && mv py/magfieldgen/.magfieldgen.py py/magfieldgen/magfieldgen.py || $(TRUE))

uninstall:
	$(RM) $(INSTALL_DIR)/$(TARGETLIB)
	cat python_install.txt | xargs rm -vrf
	$(RM) python_install.txt

rmbuild: 
	$(RM) build/$(TARGETLIB)
	($(RM) -R build || $(ECHO) "Could not remove 'build/' directory, manually remove it")
