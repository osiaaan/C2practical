# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module createdoc.sh updatePractical.sh

SUBDIRS = src doc 

doc: FORCE
	sh ./createdoc.sh
FORCE:

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune=$(DUNEROOT) CXX="$(CXX)" CC="$(CC)"

dist-hook:
	 sh ./makePatches.sh $(distdir)

