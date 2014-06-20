SHELL = /bin/sh
TOPDIR=${UQTK_SRC}
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = infer.x
SRCS = infer.cpp
OBJS = $(SRCS:.cpp=.o)
INCS = XMLreader.h posterior.h
LIBINCDIR = $(TOPDIR)/src_cpp/include
LIBBINDIR = $(TOPDIR)/src_cpp/lib
INCDIRS = -I. -I$(LIBINCDIR)
LIBS = -L$(LIBBINDIR) -luqtkmcmc -llbfgs -lUQTk -lquad -luqtktools -lcvode-2.6.0 -ldsfmt -l$(LAPACK) -l$(SLATEC) -l$(BLAS) -lxmlutils -lexpat $(FCLIB)
LDEP = $(LIBBINDIR)/libUQTk.a $(LIBBINDIR)/libquad.a $(LIBBINDIR)/libuqtktools.a $(LIBBINDIR)/libcvode-2.6.0.a \
       $(LIBBINDIR)/libdsfmt.a $(LIBBINDIR)/lib$(LAPACK).a $(LIBBINDIR)/lib$(SLATEC).a $(LIBBINDIR)/lib$(BLAS).a
DEFS = -D__$(FTNNAME)
.SUFFIXES : .cpp

all: $(TARGET)

$(TARGET): $(OBJS) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
	rm -f *.pdf *.log *.pyc *.xml
	rm -f chain.dat inputdata.dat

.f.o:
	$(F77) $(FFLAGS) $(DEFS) $(INCDIRS) -c $*.f

.c.o:
	$(CC) $(CCFLAGS) $(DEFS) $(INCDIRS) -c $*.c

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

examples: $(TARGET)
	@echo "Examples Line Inference"
	./bayesianInference1D.py

test: $(TARGET)
	@echo "The test target has not been implemented for $(TARGET)"
