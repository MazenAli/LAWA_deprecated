include $(LAWA_HOME)/Makefile.common

SUBDIRS := math aux 
SOURCES  := $(wildcard *.cc)
OBJECTS  := $(SOURCES:.cc=.o)
DEPS     := $(SOURCES:.cc=.d)
CPPFLAGS += -I $(LAWA_HOME)

DYLIB    = liblawa.$(DYLIB_EXT)

all: $(DYLIB)

$(DYLIB): $(OBJECTS)
	$(SILENT) for i in $(SUBDIRS); \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i all; \
	done;
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) $(LDFLAGS) -o $(DYLIB) *.o aux/*.o math/*.o $(LDLIBS)
	$(CP) $(DYLIB) $(LAWA_HOME)

clean:
	$(RM) *.o *.d $(DYLIB)
	$(SILENT) for i in $(SUBDIRS) ; \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i clean; \
	done;
	$(RM) -f *.tmp

# automatic dependencies -------------------------------------------------------
%.d: %.cc
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
