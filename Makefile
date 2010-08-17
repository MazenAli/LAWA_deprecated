include $(LAWA_HOME)/Makefile.common
SUBDIRS := lawa extensions

all:
	$(SILENT) for i in $(SUBDIRS); \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i all; \
	done;
	
clean:
	$(SILENT) for i in $(SUBDIRS) ; \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i clean; \
	done;
	$(RM) libextensionsflens.$(DYLIB_EXT) liblawa.$(DYLIB_EXT) \
              liblawamath.$(DYLIB_EXT)
	$(RM) *.tmp
# DO NOT DELETE
