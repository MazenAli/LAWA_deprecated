include $(LAWA_HOME)/Makefile.common

SUBDIRS := dku

SOURCEFILES  = $(wildcard *.cc)
TARGETS     = $(patsubst %.cc,%,$(SOURCEFILES))
CPPFLAGS    += -I$(FLENS_HOME) -I$(LAWA_HOME)
LDFLAGS     += -L$(LAWA_HOME) -llawamath -llawa -lextensionsgmpmath

all: $(TARGETS)
	$(SILENT) for i in $(SUBDIRS); \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i all; \
	done;
clean:
	$(RM) *.o
	$(RM) $(TARGETS)
	$(RM) *.gps *.dat *.pdf *.m *~ *.tex *.txt *.out *.png *.eps
	$(RM) -r *.dSYM
	$(SILENT) for i in $(SUBDIRS) ; \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i clean; \
	done;
	