AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

# program names 
check_PROGRAMS = laplace 
nodist_laplace_SOURCES = laplace.cc 

CLEANFILES = manager.*.log *.dgf.* dgfparser.log *.vtu corner.git*

practical7.pdf: practical7.tex
	pdflatex practical7.tex

aufgabe7dir=$(includedir)/src/7/

EXTRA_DIST= 
