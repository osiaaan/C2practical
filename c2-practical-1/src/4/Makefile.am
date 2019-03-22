AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDUNE_DEVEL_MODE
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)

# program names, create new name for spmatrix 
check_PROGRAMS = interpolate l2projection
nodist_interpolate_SOURCES = interpolate.cc  
nodist_l2projection_SOURCES = l2projection.cc  

CLEANFILES = *.log *.vtu *.aux *.dvi

practical4.pdf: practical4.tex
	pdflatex practical4.tex

aufgabe4dir=$(includedir)/src/4/

EXTRA_DIST=
