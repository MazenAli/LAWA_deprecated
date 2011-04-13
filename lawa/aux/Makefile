include $(LAWA_HOME)/Makefile.common

SOURCES  := $(wildcard *.cpp)
OBJECTS  := $(SOURCES:.cpp=.o)
DEPS     := $(SOURCES:.cpp=.d)
CPPFLAGS += -I $(LAWA_HOME) -I $(FLENS_HOME)

all: $(OBJECTS)
    
clean:
	$(RM) *.o *.d $(DYLIB)

# automatic dependencies -------------------------------------------------------
%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
