TARGET = main
CPP = g++

# 
MAIN = /ssdwork/psp_builds/rosetta_src_2018.44.60488_bundle/main
LINUXVER = 4.19
GCCVER = 8.2

RUNNAME=run.sh
# libbasic  libObjexxFCL libprotocols libcore etc
RUNLIBS1=$(MAIN)/source/build/src/release/linux/$(LINUXVER)/64/x86/gcc/$(GCCVER)/default
# libcifparse.so  libcppdb.so  libsqlite3.so  libxml2.so  libzmq.so
RUNLIBS2=$(MAIN)/source/build/external/release/linux/$(LINUXVER)/64/x86/gcc/$(GCCVER)/default
# rosetta options
RUNFLAGS=-mute all

##

CPPFLAGS = -c -std=c++11 -Wall -Wextra -Wpedantic -ffor-scope -MD
CPPFLAGSEXTRA = -pipe -pedantic -Wno-long-long -Wno-strict-aliasing -march=core2 -mtune=generic -O3 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s -Wno-unused-variable -Wno-unused-parameter -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS -DPTR_STD -DNDEBUG 

IS = -isystem $(MAIN)/source/external/boost_1_55_0/ -isystem $(MAIN)/source/external/ -isystem $(MAIN)/source/external/include/ -isystem $(MAIN)/source/external/dbio/

INCLUDE = -I$(MAIN)/source/src -I$(MAIN)/source/external/include -I$(MAIN)/source/src/platform/linux/64/gcc/$(GCCVER) -I$(MAIN)/source/src/platform/linux/64/gcc -I$(MAIN)/source/src/platform/linux/64 -I$(MAIN)/source/src/platform/linux -I$(MAIN)/source/external/boost_1_55_0 -I$(MAIN)/source/external/libxml2/include -I$(MAIN)/source/external -I$(MAIN)/source/external/dbio -I/usr/include -I/usr/local/include

LIBS1 = -L$(MAIN)/source/external/lib -L$(MAIN)/source/build/src/release/linux/$(LINUXVER)/64/x86/gcc/$(GCCVER)/default -L$(MAIN)/source/build/src/release/linux/$(LINUXVER)/64/x86/gcc/$(GCCVER)/default -L$(MAIN)/source/src -L$(MAIN)/source/build/external/release/linux/$(LINUXVER)/64/x86/gcc/$(GCCVER)/default -L$(MAIN)/external
LIBS2 = -L/usr/lib -L/usr/local/lib
LIBS3 = -ldevel -lprotocols.8 -lprotocols.7 -lprotocols_e.6 -lprotocols_d.6 -lprotocols_c.6 -lprotocols_b.6 -lprotocols_a.6 -lprotocols_h.5 -lprotocols_g.5 -lprotocols_f.5 -lprotocols_e.5 -lprotocols_d.5 -lprotocols_c.5 -lprotocols_b.5 -lprotocols_a.5 -lprotocols.4 -lprotocols.3 -lprotocols_b.2 -lprotocols_a.2 -lprotocols.1 -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3 -lcifparse -lxml2 -lzmq

SRCPATH = ./src

OBJDIR_RELEASE = obj/Release

OBJ_RELEASE = \
	$(OBJDIR_RELEASE)/main.o \
	$(OBJDIR_RELEASE)/dunbrackdata.o \
	$(OBJDIR_RELEASE)/bbutils.o \
	$(OBJDIR_RELEASE)/bbdep_sm.o 

HEADERS = \
	$(SRCPATH)/dunbrackdata.hh \
	$(SRCPATH)/bbutils.hh \
	$(SRCPATH)/bbdep_sm.hh \
	$(SRCPATH)/data_writing.hh 

all: release

clean: clean_release

release: before_release out_release

before_release:
	$(file >$(RUNNAME),export LD_LIBRARY_PATH=$(RUNLIBS1):$(RUNLIBS2))
	$(file >>$(RUNNAME),./main -database $(MAIN)/database $(RUNFLAGS))
	chmod +x $(RUNNAME)
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)
	test -d maps || mkdir -p maps

out_release: before_release $(OBJ_RELEASE) $(HEADERS)
	$(CPP) -o $(TARGET) $(LIBS1) $(OBJDIR_RELEASE)/*.o $(LIBS2) $(LIBS3)

clean_release:
	rm $(OBJDIR_RELEASE)/*.o
	rm $(OBJDIR_RELEASE)/*.d
	rm $(TARGET)
	rm $(RUNNAME)

$(OBJDIR_RELEASE)/main.o: $(SRCPATH)/main.cc
	$(CPP) $(CPPFLAGS) $(IS) $(CPPFLAGSEXTRA) $(INCLUDE) $(SRCPATH)/main.cc -o $(OBJDIR_RELEASE)/main.o

$(OBJDIR_RELEASE)/dunbrackdata.o: $(SRCPATH)/dunbrackdata.cc
	$(CPP) $(CPPFLAGS) $(IS) $(CPPFLAGSEXTRA) $(INCLUDE) $(SRCPATH)/dunbrackdata.cc -o $(OBJDIR_RELEASE)/dunbrackdata.o

$(OBJDIR_RELEASE)/bbutils.o: $(SRCPATH)/bbutils.cc
	$(CPP) $(CPPFLAGS) $(SRCPATH)/bbutils.cc -o $(OBJDIR_RELEASE)/bbutils.o

$(OBJDIR_RELEASE)/bbdep_sm.o: $(SRCPATH)/bbdep_sm.cc
	$(CPP) $(CPPFLAGS) $(IS) $(CPPFLAGSEXTRA) $(INCLUDE) $(SRCPATH)/bbdep_sm.cc -o $(OBJDIR_RELEASE)/bbdep_sm.o

-include $(OBJ_RELEASE:.o=.d)
