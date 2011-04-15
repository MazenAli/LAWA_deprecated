FLENS_HOME=/Users/stip/FLENS-lite
LAWA_HOME=/Users/stip/Desktop/LAWA-lite
include $(LAWA_HOME)/config


INCDIRS += -I$(FLENS_HOME) -I$(LAWA_HOME)

LDFLAGS += -L$(FLENS_HOME) -L$(LAWA_HOME)
LDFLAGS += -llawa -lextensionsflens -lextensionssparsegrid

EX_LDFLAGS += -bind_at_load

ifdef GMP_HOME
  CXXFLAGS += -DMPFR -I$(MPFR_HOME)/include
  LDFLAGS += -L$(MPFR_HOME)/lib -lgmpxx -lgmp
endif


SOURCEFILES  = $(wildcard *.cpp)
TARGETS      = $(patsubst %.cpp,%,$(SOURCEFILES))
LINK         = $(LDFLAGS) $(EX_LDFLAGS)

all: $(TARGETS) 

clean:
	$(RM) *.o
	$(RM) $(TARGETS)
    
$(TARGETS): % : %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCDIRS) -I . -o $@ $+ $(LINK)
	
