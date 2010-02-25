#-------------------------------------------------------------------------
#  File         :  Makefile.in
#-------------------------------------------------------------------------
#
#  Purpose      :  Compilation
#
#  Remarks      :
#
#  References   :
#
#  Revisions    :
#-------------------------------------------------------------------------
#  Parallel Particle Mesh Library (PPM)
#  Institute of Computational Science
#  ETH Zentrum, Hirschengraben 84
#  CH-8092 Zurich, Switzerland
#-------------------------------------------------------------------------

# These variables are set when configure is ran
prefix = /Users/omer/ppm/ppm/libppmnumerics
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
builddir = .
LIBS = -lmpi -lppm -lm 
LDFLAGS = -L/Users/omer/metis-4.0/lib -L../libppm/lib
CFLAGS = -g -O2
FCLIBS =  -L/Users/omer/metis-4.0/lib -L/opt/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -lutil -L/usr/bin/ifort-11.1-base/lib -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64/ -L/usr/lib/i686-apple-darwin10/4.2.1/ -L/usr/lib/ -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/ -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/../../../i686-apple-darwin10/4.2.1/ -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/../../.. /usr/bin/ifort-11.1-base/lib/libifport.a /usr/bin/ifort-11.1-base/lib/libifcore.a /usr/bin/ifort-11.1-base/lib/libimf.a /usr/bin/ifort-11.1-base/lib/libsvml.a /usr/bin/ifort-11.1-base/lib/libipgo.a -lSystemStubs -lmx /usr/bin/ifort-11.1-base/lib/libirc.a -lpthread -ldl
FCFLAGS = -g -FR -I../libppm/include
FC = mpif90
CC = mpicc
CXX = mpic++
DEBUG = 

# These variables are standard
LIB_PPM := libppmnumberics.a
TARGET  := $(builddir)/lib/$(LIB_PPM)
SHELL := /bin/sh
CPP := cpp

# These are important build directories
SRC_DIR := $(builddir)/src
OBJ_DIR := $(builddir)/objects
MODULES_DIR := $(builddir)/include

# These are paths that get included during compilation
VPATH := $(SRC_DIR):
VPATH += $(MODULES_DIR):
INCLS := $(patsubst %,-I%,$(subst :, ,$(VPATH)))

# These are the files that get generated and used during compilation
SOURCES := $(notdir $(wildcard $(SRC_DIR)/ppm_module_*.f))
OBJECTS := $(SOURCES:%.f=$(OBJ_DIR)/%.o)
MODULES := $(SOURCES:%.f=$(MODULES_DIR)/%.mod)
MODSRCS := $(SOURCES:%.f=$(MODULES_DIR)/__%.f)
DEPENDENCIES := $(SOURCES:%.f=$(OBJ_DIR)/%.d)

# This creates the install directories if they don't exist
$(warning Checking for directories...)
$(shell mkdir $(OBJ_DIR))
$(shell mkdir $(MODULES_DIR))
$(shell mkdir $(builddir)/lib)
$(shell mkdir $(libdir))

.DEFAULT: ;

all: $(TARGET)

# This archives all of the objects in the PPM library
$(TARGET): $(OBJECTS)
	ar crus $@ $(OBJECTS)

# This creates the file dependencies
# 1) we use the preprocessor to find the includes
# 2) we add the dependency file as a target
# 3) find INCLUDE and USE statements that are not inside a comment
$(OBJ_DIR)/%.d: $(SRC_DIR)/%.f
	@echo '[ quietly making' $@ ']'
	@$(CPP) $(INCLS) -M $< | \
	sed -e 's#$*.o#$(OBJ_DIR)/$*.o $(OBJ_DIR)/$*.d#' \
		-e 's#$$#\\#' -e's#\\\\#\\#' > $@
	@$(CPP) -P $(INCLS) $< > $(OBJ_DIR)/__$*.f
	@grep "INCLUDE " $(OBJ_DIR)/__$*.f | \
		sed -e 's#^[ \t]*##;s#[ \t]*##' \
			-e  '/^!/d;/mpif.h/d' \
			-e  '/^!/d;/fftw3.f/d' \
			-e 's#INCLUDE ##;s#$$#\\#' \
			-e 's#"##g' -e "s#'##g" -e 's# ##g' >> $@
	@echo '# end of source dependencies for .o and .d files' >> $@
	@echo ''$(OBJ_DIR)/$*.o ':\' >> $@ 
	@grep "USE " $(OBJ_DIR)/__$*.f | \
		sed -e 's#^[ \t]*##;s#[ \t]*##' \
			-e '/^!/d' \
			-e 's#,.*##' | \
		sed -e 's#USE #$(OBJ_DIR)/#' \
			-e 's# ##g;s#$$#.o\\#' >> $@
	@echo '# end of module dependencies for .o file' >> $@
	@rm $(OBJ_DIR)/__$*.f

# This handles the pre-processing and does the actual compiling
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(CPP) -P $(INCLS) $< > $(OBJ_DIR)/__$*.f
	$(FC) $(INCLS) $(LDFLAGS) $(FCFLAGS) $(LIBS) $(DEBUG) -Llib -c -o $@ $(OBJ_DIR)/__$*.f
	@mv $(builddir)/*.mod $(MODULES_DIR)
	@rm $(OBJ_DIR)/__$*.f

# This is used to clean up the files created when running make
clean:
	rm -fR $(OBJ_DIR)
	rm -fR $(MODULES_DIR)
	rm -fR $(builddir)/lib
	rm -f $(libdir)/$(LIB_PPM)

# This copies the PPM library into libdir
install: all
	@echo '[ deploying to '$(libdir)']'
	@cp $(TARGET) $(libdir)
	@ranlib $(libdir)/$(LIB_PPM)

# This cleans, compiles, and copies the PPM library
new: clean all install

# This ensures all dependency files are up-to-date
include $(DEPENDENCIES)
