include $(LAWA_HOME)/Makefile.common

DYLIB    = libextensionsflens.$(DYLIB_EXT)

SOURCES  := $(wildcard *.cpp)
OBJECTS  := $(SOURCES:.cpp=.o)
DEPS     := $(SOURCES:.cpp=.d)
CPPFLAGS += -I $(LAWA_HOME)

all: $(DYLIB)

$(DYLIB): $(OBJECTS)
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) $(LDFLAGS) -o $(DYLIB) *.o $(LDLIBS)
	$(CP) $(DYLIB) $(LAWA_HOME)

clean:
	$(RM) *.o *.d $(DYLIB)

# automatic dependencies -------------------------------------------------------
%.d: %.cpp
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
