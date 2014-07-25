SHELL = /bin/sh
TOPDIR=${UQTK_SRC}
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = infer.x
SRCS1 = infer.cpp 
OBJS1 = $(SRCS1:.cpp=.o)
SRCS2 = infer2D.cpp 
OBJS2 = $(SRCS2:.cpp=.o)
INCS = XMLreader.h posterior.h XMLreader2D.h posterior2D.h
LIBINCDIR = $(TOPDIR)/src_cpp/include
LIBBINDIR = $(TOPDIR)/src_cpp/lib
INCDIRS = -I. -I$(LIBINCDIR)
LIBS = -L$(LIBBINDIR) -luqtkmcmc -llbfgs -lUQTk -lquad -luqtktools -lcvode-2.6.0 -ldsfmt -l$(LAPACK) -l$(SLATEC) -l$(BLAS) -lxmlutils -lexpat $(FCLIB)
LDEP = $(LIBBINDIR)/libUQTk.a $(LIBBINDIR)/libquad.a $(LIBBINDIR)/libuqtktools.a $(LIBBINDIR)/libcvode-2.6.0.a \
       $(LIBBINDIR)/libdsfmt.a $(LIBBINDIR)/lib$(LAPACK).a $(LIBBINDIR)/lib$(SLATEC).a $(LIBBINDIR)/lib$(BLAS).a
DEFS = -D__$(FTNNAME)
.SUFFIXES : .cpp

all: infer.x infer2D.x

infer.x: $(OBJS1) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o infer.x $(OBJS1) $(LIBS)

infer2D.x: $(OBJS2) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o infer2D.x $(OBJS2) $(LIBS)

clean:
	rm -f $(OBJS1) $(OBJS2) infer.x infer2D.x
	rm -f *.pdf *.log *.pyc *.xml
	rm -f chain.dat inputdata.dat

.f.o:
	$(F77) $(FFLAGS) $(DEFS) $(INCDIRS) -c $*.f

.c.o:
	$(CC) $(CCFLAGS) $(DEFS) $(INCDIRS) -c $*.c

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

test: $(TARGET)
	@echo "The test target has not been implemented for $(TARGET)"
