#-------------------------------------------------------------------------
#  File         :  Makefile.in
#-------------------------------------------------------------------------
#
#  Purpose      :  Makefile template
#
#  Remarks      :
#
#  References   :
#
#  Revisions    :
#-------------------------------------------------------------------------
#  Parallel Particle Mesh Library (PPM)
#  ETH Zurich
#  CH-8092 Zurich, Switzerland
#-------------------------------------------------------------------------

# These variables are set when configure is ran
prefix = /usr/local/ppm/build
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
builddir = .
LIBS = -lppm -lppm -lm 
LDFLAGS = -L/usr/local/metis-4.0_gcc -L/home/jtra/phd/naga/src/libppm/lib -L/home/jtra/phd/naga/src/libppm/lib
CFLAGS = -O3
FCLIBS =  -L/usr/local/metis-4.0_gcc -L/usr/local/openmpi-1.4.3_gcc/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3/../../.. -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -lgfortranbegin -lgfortran -lm -lpthread
FCFLAGS = -O3 -ffree-form -I/home/jtra/phd/naga/src/libppm/include
FC = mpif90
CC = mpicc
CXX = mpic++
DEBUG = 
MODFLAG = -I

# These variables are standard
LIB_PPM := libppmnumerics.a
TARGET  := $(builddir)/lib/$(LIB_PPM)
SHELL := /bin/sh
CPP := cpp

# These are important build directories
SRC_DIR := $(builddir)/src
OBJ_DIR := $(builddir)/objects
MODULES_DIR := $(builddir)/include

# These are paths that get included during pre-processing
CPPVPATH := $(SRC_DIR):
CPPVPATH += $(MODULES_DIR):
CPPINCLS := $(patsubst %,-I%, $(subst :, ,$(CPPVPATH)))

# These are paths that get included during compilation
VPATH := $(patsubst %,-I%, $(subst :, ,$(SRC_DIR))):
VPATH += $(patsubst %,$(MODFLAG)%, $(subst :, ,$(MODULES_DIR))):
INCLS := $(subst :, ,$(VPATH))

# These are the files that get generated and used during compilation
SOURCES := $(notdir $(wildcard $(SRC_DIR)/ppm_module_*.f))
OBJECTS := $(SOURCES:%.f=$(OBJ_DIR)/%.o)
MODULES := $(SOURCES:%.f=$(MODULES_DIR)/%.mod)
MODSRCS := $(SOURCES:%.f=$(MODULES_DIR)/__%.f)
DEPENDENCIES := $(SOURCES:%.f=$(OBJ_DIR)/%.d)

# This creates the install directories if they don't exist
$(warning Checking for directories...)
$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))
$(shell test -d $(MODULES_DIR) || mkdir $(MODULES_DIR))
$(shell test -d $(builddir)/lib || mkdir $(builddir)/lib)
$(shell test -d $(libdir) || mkdir $(libdir))
$(warning done.)

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
	@$(CPP) $(CPPINCLS) -M $< | \
        sed -e 's#$*.o#$(OBJ_DIR)/$*.o $@#' \
            -e 's#$$#\\#' \
            -e 's#\\\\#\\#' > $@
	@$(CPP) -P $(CPPINCLS) $< > $(OBJ_DIR)/__$*.f
	@grep "INCLUDE " $(OBJ_DIR)/__$*.f | \
        sed -e 's#^[ \t]*##' \
            -e '/^!/d' \
            -e '/mpif.h/d' \
            -e 's#INCLUDE #$(SRC_DIR)/#' \
            -e 's#$$# \\#' \
            -e "s#'##g" >> $@
	@echo '# end of source dependencies for .o and .d files' >> $@
	@echo ''$(OBJ_DIR)/$*.o': \' >> $@
	@grep "USE " $(OBJ_DIR)/__$*.f | \
        sed -e 's#^[ \t]*##' \
            -e '/^!/d' \
            -e 's#,.*##' \
            -e 's#USE #$(OBJ_DIR)/#' \
            -e 's#$$#.o \\#' >> $@
	@echo '# end of module dependencies for .o file' >> $@
	@rm $(OBJ_DIR)/__$*.f

# This handles the pre-processing and does the actual compiling
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(CPP) -P $(CPPINCLS) $< > $(OBJ_DIR)/__$*.f
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
# WARNING: The directive below will fail silently, in case of problems
# remove the leading s: include $(DEPENDENCIES)
sinclude $(DEPENDENCIES)
