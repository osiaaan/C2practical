AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

check_PROGRAMS = laplace  l2projection 
nodist_laplace_SOURCES = laplace.cc 
nodist_l2projection_SOURCES = l2projection.cc 

CLEANFILES = manager.*.log *.dgf.* dgfparser.log *.vtu corner.git*

aufgabe8dir=$(includedir)/src/8/

EXTRA_DIST= 
