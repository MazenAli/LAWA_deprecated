include $(LAWA_HOME)/config


INCDIRS += -I$(FLENS_HOME) -I$(LAWA_HOME)
LDFLAGS += -L$(FLENS_HOME) -L$(LAWA_HOME)
EX_LDFLAGS += -bind_at_load


SOURCEFILES  = $(wildcard *.cpp)
TARGETS      = $(patsubst %.cpp,%,$(SOURCEFILES))
LINK         = $(LDFLAGS) $(EX_LDFLAGS)

all: $(TARGETS) 

clean:
	$(RM) *.o
	$(RM) $(TARGETS)
    
$(TARGETS): % : %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -I . -o $@ $+ $(LINK)
	
