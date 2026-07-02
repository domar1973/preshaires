FC ?= gfortran
FFLAGS ?= -O2 -Wall

ROOT := ../..
NR := $(ROOT)/src/NR
PROBE := dbskr3_probe
SOURCES := dbskr3_probe.f $(NR)/dbskr3_fromNR.f $(NR)/bessik.f $(NR)/beschb.f $(NR)/chebev.f

.PHONY: all run clean

all: $(PROBE)

$(PROBE): $(SOURCES)
	$(FC) $(FFLAGS) -o $@ $(SOURCES)

run: $(PROBE)
	./$(PROBE)

clean:
	rm -f $(PROBE)
