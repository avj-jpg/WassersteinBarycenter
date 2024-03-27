-include Make.user

ifndef MFEM_DIR
	$(error MFEM_DIR is not set. Please create a Make.user file to set MFEM_DIR)
endif

ifndef MFEM_BUILD_DIR
	MFEM_BUILD_DIR=$(MFEM_DIR)
endif

CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk

-include $(CONFIG_MK)

# Overall structure
KRON_BUILD_DIR=build/kron
KRON_SOURCES=kron_mult.cpp kron_quad.cpp
KRON_OBJECTS=$(KRON_SOURCES:%.cpp=$(KRON_BUILD_DIR)/%.o)

BUILD_DIR=build
SOURCES= bary_rhs.cpp  wass_laplace.cpp wass_multigrid.cpp 
OBJECTS=$(SOURCES:%.cpp=$(BUILD_DIR)/%.o)

APPS= bary_pdhg surf_pdhg 
APP_SRC=$(APPS:%=%.cpp)

CXXFLAGS=-g $(MFEM_CXXFLAGS)

# Use compiler configuration from MFEM
LFLGAS=$(MFEM_LIBS)
INCFLAGS=$(MFEM_INCFLAGS) -I kron

.PHONY: all clean style deps infodir
all: $(APPS)

# Build the executable
$(APPS):%: $(BUILD_DIR)/%.o $(OBJECTS) $(KRON_OBJECTS) $(MFEM_LIB_FILE)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) $< $(OBJECTS) $(KRON_OBJECTS) $(LFLGAS) -o $@

$(BUILD_DIR)/%.o: %.cpp makefile Make.user | $(BUILD_DIR)
	$(MFEM_CXX) -c $(CXXFLAGS) $(INCFLAGS) -o $@ $<

$(KRON_BUILD_DIR)/%.o: kron/%.cpp makefile Make.user | $(KRON_BUILD_DIR)
	$(MFEM_CXX) -c $(CXXFLAGS) $(INCFLAGS) -o $@ $<

# Ensure build directory exists
$(BUILD_DIR):
	mkdir -p $@

# Ensure kron directory exists
$(KRON_BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(EXE) $(BUILD_DIR) $(KRON_BUILD_DIR)

infodir:
	$(info $(MFEM_DIR))
	@true

FORMAT_FILES = $(wildcard *.?pp)
ASTYLE = astyle --options=$(MFEM_DIR)/config/mfem.astylerc
style:
	@if ! $(ASTYLE) $(FORMAT_FILES) | grep Formatted; then\
	   echo "No source files were changed.";\
	fi
