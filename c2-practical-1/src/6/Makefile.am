AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

# program names 
check_PROGRAMS = l2projection laplace
nodist_l2projection_SOURCES = l2projection.cc 
nodist_laplace_SOURCES = laplace.cc 

CLEANFILES = manager.*.log *.dgf.* dgfparser.log *.vtu

practical6.pdf: practical6.tex
	pdflatex practical6.tex

aufgabe6dir=$(includedir)/src/6/

EXTRA_DIST= 
