TARGET = main
CPP = g++
MAIN = /ssdwork/psp_builds/rosetta_src_2018.42.60459_bundle/main

LINUXVER = 4.18
GCCVER = 8.2

RUNNAME=run.sh
# libbasic  libObjexxFCL libprotocols libcore etc
RUNLIBS1=$(MAIN)/source/build/src/release/linux/$(LINUXVER)/64/x86/gcc/$(GCCVER)/default
# libcifparse.so  libcppdb.so  libsqlite3.so  libxml2.so  libzmq.so
RUNLIBS2=$(MAIN)/source/build/external/release/linux/$(LINUXVER)/64/x86/gcc/$(GCCVER)/default
# rosetta options
RUNFLAGS=-mute all

all:
	$(file >$(RUNNAME),export LD_LIBRARY_PATH=$(RUNLIBS1):$(RUNLIBS2))
	$(file >>$(RUNNAME),./main -database $(MAIN)/database $(RUNFLAGS))