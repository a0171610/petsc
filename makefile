CC=gcc
FC=gfortran
CFLAGS=-Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fno-stack-check -fvisibility=hidden -g3
FFLAGS=-Wl,-bind_at_load -Wl,-multiply_defined,suppress -Wl,-multiply_defined -Wl,suppress -Wl,-commons,use_dylibs -Wl,-search_paths_first -Wl,-no_compact_unwind  -fPIC -Wall -ffree-line-length-0 -Wno-lto-type-mismatch -Wno-unused-dummy-argument -g -O0   -fPIC -Wall -ffree-line-length-0 -Wno-lto-type-mismatch -Wno-unused-dummy-argument -g -O0
#this 'CPPFLAGS' is include option 
#CPPFLAGS=-I/Users/okazakikouhei/petsc/include -I/Users/okazakikouhei/petsc/arch-darwin-c-debug/include
#LDFLAGS=-Wl,-rpath,/Users/okazakikouhei/petsc/arch-darwin-c-debug/lib -L/Users/okazakikouhei/petsc/arch-darwin-c-debug/lib -Wl,-rpath,/Users/okazakikouhei/petsc/arch-darwin-c-debug/lib -L/Users/okazakikouhei/petsc/arch-darwin-c-debug/lib
CPPFLAGS=-I/usr/local/Cellar/petsc/3.17.2/include
LDFLAGS=-L/usr/local/Cellar/petsc/3.17.2/lib
LDLIBS=-lpetsc -lm

.SUFFIXES :
.SUFFIXES : .F90 .o

SRC = sparse_module.F90 main.F90
OBJ = ${SRC:.F90=.o}
TARGET=adv

$(TARGET) : $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LDLIBS) -o $@

main.o: sparse_module.o
solver_module.o : user_module.o
example.o : user_module.o solver_module.o

SRC_EXAMPLE = solver_module.F90 user_module.F90 example.F90
OBJ_EXAMPLE = ${SRC_EXAMPLE:.F90=.o}
TARGET_EXAMPLE=sparse

$(TARGET_EXAMPLE) : $(OBJ_EXAMPLE)
	$(FC) $(LDFLAGS) $(OBJ_EXAMPLE) $(LDLIBS) -o $@

clean :
	rm -f *.o *.mod $(TARGET) *.dat $(TARGET).log $(TARGET_EXAMPLE)

.F90.o :
	$(FC) $(CPPFLAGS) $< -c