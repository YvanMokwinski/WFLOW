DISTRIB=$(PWD)
include $(DISTRIB)/config.mk
OPTIMIZE_FLAGS=-O2 -DNDEBUG #-Wall -Werror
DEBUG_FLAGS=-g #-Wall -Werror
all: 


# OPTIMIZED VERSION
	mkdir -p $(PLATFORM)_$(CC)
	+make $@ -C $(PLATFORM)_$(CC)  -f $(DISTRIB)/src/Makefile DISTRIB=$(DISTRIB) DEBUG="" FLAGS="$(OPTIMIZE_FLAGS)"


# DEBUG VERSION
	mkdir -p $(PLATFORM)_$(CC)Debug
	+make $@ -C $(PLATFORM)_$(CC)Debug -f $(DISTRIB)/src/Makefile DISTRIB=$(DISTRIB) DEBUG=Debug FLAGS="$(DEBUG_FLAGS)"


help:
	@echo ""
	@echo "Makefile for MNS"
	@echo "written by Yvan Mokwinski"
	@echo ""
	@echo ""
	@echo "make all"
	@echo "     build optimized and debug programs"
	@echo ""
	@echo "make build_doc"
	@echo "     build code documentation using Doxygen"
	@echo "     the  main page of the generated documentation is "
	@echo "     doc/html/html/index.html"
	@echo ""
	@echo "make build_report"
	@echo "     build the report for the practicum"
	@echo ""
	@echo "make clean"
	@echo "     clean object files"
	@echo ""
	@echo "make realclean"
	@echo "     clean object files and Doxygen documentation"
	@echo ""
	@echo ""

build_doc:
	cd doc;doxygen Doxyfile

clean_doc:
	rm -rf doc/html

build_report:
	cd report;make $@

clean_report:
	cd report;make $@

cleandistrib:
	rm -rf $(PLATFORM)_$(CC)
	rm -rf $(PLATFORM)_$(CC)Debug
	rm -f  lib/*.a
	rm -f  lib/*~ *~ src/*~

clean:cleandistrib

realclean:cleandistrib clean_doc clean_report


