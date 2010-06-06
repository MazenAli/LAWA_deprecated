include $(LAWA_HOME)/Makefile.common

SOURCEFILES  = $(wildcard *.cc)
TARGETS     = $(patsubst %.cc,%,$(SOURCEFILES))
CPPFLAGS    += -I$(FLENS_HOME) -I$(LAWA_HOME)
LDFLAGS     += -L$(LAWA_HOME) -llawamath -llawa

all: $(TARGETS)

clean:
	$(RM) *.o
	$(RM) $(TARGETS)
	$(RM) *.gps *.dat *.pdf *.m *~ *.tex *.txt *.out *.png *.eps
	$(RM) -r *.dSYM
