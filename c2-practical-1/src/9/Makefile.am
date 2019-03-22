AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

check_PROGRAMS = adapt surfell surfparab elastic3D elastic2D laplace3D
nodist_adapt_SOURCES = adapt.cc
nodist_surfell_SOURCES = surfell.cc 
nodist_surfparab_SOURCES = surfparab.cc
nodist_elastic2D_SOURCES = elastic2D.cc
nodist_elastic3D_SOURCES = elastic3D.cc
nodist_laplace3D_SOURCES = laplace3D.cc

EXTRA_DIST= 
CLEANFILES = manager.*.log *.dgf.* dgfparser.log *.vtu corner.git*

aufgabe9dir=$(includedir)/src/9/


