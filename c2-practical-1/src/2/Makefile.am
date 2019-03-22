AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

check_PROGRAMS = exercise1d exercise2
nodist_exercise2_SOURCES = exercise2.cc
nodist_exercise1d_SOURCES = exercise1d.cc

CLEANFILES = manager.*.log *.dgf.* dgfparser.log *.vtu *.gnu

practical2.pdf: practical2.tex
	pdflatex practical2.tex

aufgabe2dir=$(includedir)/src/2/

EXTRA_DIST=
