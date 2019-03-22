AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

# program names 
check_PROGRAMS = l2projection
nodist_l2projection_SOURCES = l2projection.cc 

CLEANFILES = manager.*.log *.dgf.* dgfparser.log *.vtu

practical5.pdf: practical5.tex
	pdflatex practical5.tex

aufgabe5dir=$(includedir)/src/5/

EXTRA_DIST=
