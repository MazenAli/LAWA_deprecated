include $(LAWA_HOME)/Makefile.common

SOURCES  := $(wildcard *.cc)
OBJECTS  := $(SOURCES:.cc=.o)
DEPS     := $(SOURCES:.cc=.d)
CPPFLAGS += -I $(LAWA_HOME) -I $(FLENS_HOME)

all:
	
clean:
	$(RM) *.o *.d $(DYLIB)

# automatic dependencies -------------------------------------------------------
%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
