include $(LAWA_HOME)/Makefile.common

DYLIB    = libextensionssparsegrid.$(DYLIB_EXT)

SOURCES  := sandia_rules.cc sparse_grid_mixed.cc sparse_grid_open.cc
OBJECTS  := $(SOURCES:.cc=.o)
DEPS     := $(SOURCES:.cc=.d)

all: $(DYLIB)

$(DYLIB): $(OBJECTS)
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) $(LDFLAGS) -o $(DYLIB) $(OBJECTS) $(LDLIBS)
	$(CP) $(DYLIB) $(LAWA_HOME)

clean:
	$(RM) *.o *.d $(DYLIB)

# automatic dependencies -------------------------------------------------------
%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
