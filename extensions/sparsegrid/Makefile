include $(LAWA_HOME)/Makefile.common

DYLIB    = libextensionssparsegrid.$(DYLIB_EXT)

SOURCES  := $(wildcard *.C)
OBJECTS  := $(SOURCES:.C=.o)
DEPS     := $(SOURCES:.C=.d)

all: $(DYLIB)

$(DYLIB): $(OBJECTS)
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) $(LDFLAGS) -o $(DYLIB) *.o $(LDLIBS)
	$(CP) $(DYLIB) $(LAWA_HOME)

clean:
	$(RM) *.o *.d $(DYLIB)

# automatic dependencies -------------------------------------------------------
%.d: %.C
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
