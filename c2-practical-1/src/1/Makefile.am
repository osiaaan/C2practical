AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

CLEANFILES = manager.*.log *.dgf.* dgfparser.log *.vtu *.gnu

practical1.pdf: practical1.tex
	pdflatex practical1.tex

check_PROGRAMS = exercise1d
nodist_exercise1d_SOURCES = exercise1d.cc

aufgabe1dir=$(includedir)/src/1/

EXTRA_DIST=practical1.pdf *.c *.cc

