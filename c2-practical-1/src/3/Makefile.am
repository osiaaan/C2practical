AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

check_PROGRAMS = exercise2 interpolate
nodist_exercise2_SOURCES = exercise2.cc
nodist_interpolate_SOURCES = interpolate.cc 

CLEANFILES = *.vtu *.aux *.log *.dvi

practical3.pdf: practical3.tex
	pdflatex practical3.tex

aufgabe3dir=$(includedir)/src/3/
EXTRA_DIST=
