
SUBDIRS = src \
	doc \
	test_ester

EXTRA_DIST = autogen.sh

# Add an item for each and every test target to the CLEANFILES variable
CLEANFILES = *~ 

etags:
	etags `find . -name "*.[ch]*"`

ctags:
	ctags `find . -name "*.[ch]*"`

doxygen: 
	doxygen ./doc/doxygen.cfg

