# --------------------------------------------------
# Build configuration
# --------------------------------------------------

PYTHIA  ?= $(HOME)/Desktop/projects/pythia8318
FASTJET ?= $(HOME)/Desktop/projects/fastjet-install

CXX      ?= g++
ROOT     ?= root
CXXFLAGS ?= -O2 -std=c++17
CXXFLAGS += $(shell root-config --cflags) \
            -I$(PYTHIA)/include -I$(FASTJET)/include

LDFLAGS += $(shell root-config --libs) \
           -L$(PYTHIA)/lib -lpythia8 \
           -L$(FASTJET)/lib -lfastjet \
           -Wl,-rpath,$(PYTHIA)/lib \
           -Wl,-rpath,$(FASTJET)/lib

# --------------------------------------------------
# Executables
# --------------------------------------------------

PROGRAMS = jets-charm jets-inc \
           buildtree buildtree_os \
           features_jetbins features_jetbins_os \
           reconstruct-d0

HEADERS = $(wildcard macros/*.h)

.PHONY: all clean
all: $(PROGRAMS)

# All executables except reconstruct-d0 use macros/<target>.cc.
$(filter-out reconstruct-d0,$(PROGRAMS)): %: macros/%.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

reconstruct-d0: macros/reconstructD0.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# --------------------------------------------------
# Data generation
# --------------------------------------------------

EVENTS ?= 1000000
JOBS   ?= 5

# The inclusive generation script writes this file.
PARTICLE_INPUT ?= data/d0_inc_5_100.root

.PHONY: generate generate-inclusive generate-charm
generate: generate-inclusive

generate-inclusive: jets-inc
	@mkdir -p data logs
	./scripts/run_generation.sh inc $(EVENTS) $(JOBS)

generate-charm: jets-charm
	@mkdir -p data logs
	./scripts/run_generation.sh charm $(EVENTS) $(JOBS)

# Allow candidate production to generate its default input when missing.
data/d0_inc_5_100.root: | jets-inc
	@mkdir -p data logs
	./scripts/run_generation.sh inc $(EVENTS) $(JOBS)

# --------------------------------------------------
# Candidate building and feature plots
# --------------------------------------------------

BDT_INPUT  ?= data/candidates_inc_5_100.root
BDT_OUTPUT ?= tmva
BUILD_OPTS ?=

$(BDT_INPUT): $(PARTICLE_INPUT) | buildtree
	@mkdir -p data logs
	./buildtree $(PARTICLE_INPUT) --out=$@ $(BUILD_OPTS) \
	  > logs/buildtree_bdt.log 2>&1

.PHONY: candidates features-plots
candidates: $(BDT_INPUT)

FEATURE_INPUT ?= $(BDT_INPUT)
features-plots: features_jetbins $(BDT_INPUT)
	./features_jetbins $(FEATURE_INPUT)

# Optional same-sign (OS-only) candidate production and feature plots.
OS_PARTICLE_INPUT ?= $(PARTICLE_INPUT)
OS_CANDIDATE_INPUT ?= data/candidates_os_inc_5_100.root

$(OS_CANDIDATE_INPUT): $(OS_PARTICLE_INPUT) | buildtree_os
	@mkdir -p data logs
	./buildtree_os $(OS_PARTICLE_INPUT) --out=$@ \
	  > logs/buildtree_bdt_os.log 2>&1

.PHONY: candidates-os features-plots-os
candidates-os: $(OS_CANDIDATE_INPUT)

features-plots-os: features_jetbins_os $(OS_CANDIDATE_INPUT)
	./features_jetbins_os $(OS_CANDIDATE_INPUT)

# --------------------------------------------------
# TMVA training
# --------------------------------------------------

.PHONY: check-tmva train-bdt-jetbins
check-tmva:
	@root-config --features 2>/dev/null | grep -qi tmva || { \
	  echo "ERROR: this ROOT installation does not advertise TMVA support."; \
	  echo "Check that ROOT was built with TMVA and root-config is on PATH."; \
	  exit 1; \
	}

train-bdt-jetbins: check-tmva train_bdt_jetbins.C $(BDT_INPUT)
	@mkdir -p $(BDT_OUTPUT)
	$(ROOT) -l -b -q 'train_bdt_jetbins.C("$(BDT_INPUT)","$(BDT_OUTPUT)")'

# --------------------------------------------------
# Cleaning
# --------------------------------------------------

.PHONY: clean clean-bdt clean-candidates
clean:
	rm -f $(PROGRAMS)

clean-bdt:
	rm -rf $(BDT_OUTPUT) dataset_BDT_*

clean-candidates:
	rm -f $(BDT_INPUT) $(OS_CANDIDATE_INPUT)
