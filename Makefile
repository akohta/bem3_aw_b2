# makefile for GNU C compiler
CC     =gcc
CFLAGS =-O3 -Wall
LDFLAGS=-lpng -lz -lgsl -lgslcblas
MPFLAGS=-fopenmp
# link mkl lapack library, static linking 
LDF_MKL =-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl#mkl lapack link library (static linking) 
# mkl compile option
COP_MKL =-DMKL_ILP64 -m64 -I${MKLROOT}/include
SRCDIR1=d3b2_src
SRCDIR2=maw_src
SRCDIR3=com_src
OBJDIR =obj
TARGET1=aw_d3b2_create_matrix
TRGSRC1=aw_d3b2_create_matrix.c
TARGET2=aw_d3b2_bv_solver
TRGSRC2=aw_d3b2_bv_solver.c
TARGET3=example1.out
TRGSRC3=example1.c
TARGET4=example2.out
TRGSRC4=example2.c
TARGET5 =example3.out
TRGSRC5 =example3.c
TARGET6 =example4.out
TRGSRC6 =example4.c

SRCS1=$(wildcard $(SRCDIR1)/*.c)
OBJS1=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS1)) ))
HEAD1=$(wildcard $(SRCDIR1)/*.h)

SRCS2=$(wildcard $(SRCDIR2)/*.c)
OBJS2=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS2)) ))
HEAD2=$(wildcard $(SRCDIR2)/*.h)

SRCS3=$(wildcard $(SRCDIR3)/*.c)
OBJS3=$(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(notdir $(SRCS3)) ))
HEAD3=$(wildcard $(SRCDIR3)/*.h)

TRGOBJ1=$(filter-out $(OBJDIR)/$(TARGET2).o, $(OBJS1)) $(OBJS2) $(OBJS3) 
TRGOBJ2=$(filter-out $(OBJDIR)/$(TARGET1).o, $(OBJS1)) $(OBJS2) $(OBJS3) 
TRGOBJ3=$(filter-out $(OBJDIR)/$(TARGET1).o $(OBJDIR)/$(TARGET2).o, $(OBJS1)) $(OBJS2) $(OBJS3) $(patsubst %.c,%.o,$(TRGSRC3))
TRGOBJ4=$(filter-out $(OBJDIR)/$(TARGET1).o $(OBJDIR)/$(TARGET2).o, $(OBJS1)) $(OBJS2) $(OBJS3) $(patsubst %.c,%.o,$(TRGSRC4))
TRGOBJ5=$(filter-out $(OBJDIR)/$(TARGET1).o $(OBJDIR)/$(TARGET2).o, $(OBJS1)) $(OBJS2) $(OBJS3) $(patsubst %.c,%.o,$(TRGSRC5))
TRGOBJ6=$(filter-out $(OBJDIR)/$(TARGET1).o $(OBJDIR)/$(TARGET2).o, $(OBJS1)) $(OBJS2) $(OBJS3) $(patsubst %.c,%.o,$(TRGSRC6))

all : directories $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) 

directories:
	@mkdir -p $(OBJDIR)

$(TARGET1) : $(TRGOBJ1)  
	$(CC) $(CFLAGS) $^ -o $@ $(LDF_MKL) $(LDFLAGS) $(MPFLAGS)

$(TARGET2) : $(TRGOBJ2)  
	$(CC) $(CFLAGS) $^ -o $@ $(LDF_MKL) $(LDFLAGS) $(MPFLAGS)

$(TARGET3) : $(TRGOBJ3)  
	$(CC) $(CFLAGS) $^ -o $@ $(LDF_MKL) $(LDFLAGS) $(MPFLAGS)

$(TARGET4) : $(TRGOBJ4)  
	$(CC) $(CFLAGS) $^ -o $@ $(LDF_MKL) $(LDFLAGS) $(MPFLAGS)

$(TARGET5) : $(TRGOBJ5)  
	$(CC) $(CFLAGS) $^ -o $@ $(LDF_MKL) $(LDFLAGS) $(MPFLAGS)

$(TARGET6) : $(TRGOBJ6)  
	$(CC) $(CFLAGS) $^ -o $@ $(LDF_MKL) $(LDFLAGS) $(MPFLAGS)

$(OBJDIR)/%.o : $(SRCDIR1)/%.c
	$(CC) $(CFLAGS) $(COP_MKL) -I$(SRCDIR2) -I$(SRCDIR3) -c $< -o $@ $(MPFLAGS)

$(OBJDIR)/%.o : $(SRCDIR2)/%.c
	$(CC) $(CFLAGS) $(COP_MKL) -I$(SRCDIR3) -c $< -o $@ $(MPFLAGS)

$(OBJDIR)/%.o : $(SRCDIR3)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

.c.o :
	$(CC) $(CFLAGS) $(COP_MKL) -I$(SRCDIR1) -I$(SRCDIR2) -I$(SRCDIR3) -c $< $(MPFLAGS)

clean:
	@rm -rf $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(OBJDIR) ./*.o

$(OBJS1) : $(HEAD1) $(HEAD2) $(HEAD3) 
$(OBJS2) : $(HEAD2) $(HEAD3) 
$(OBJS3) : $(HEAD3) 
