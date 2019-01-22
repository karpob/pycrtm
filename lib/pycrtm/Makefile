# -*- makefile -*
#
#	${F2PY} ${F2PY_FLAGS} -c -L${LIB} -lcrtm -lgomp -I${INC} wrapped.o -m ${MODULE}  $< only: run_forward run_k_matrix
#  FORT = ifort
#  F2PY = f2py --fcompiler=intelem --f90flags='${FCFLAGS}'
  F2PY = f2py --fcompiler=${F2PY_COMPILER} --f90flags='${FCFLAGS}' #--f90exec=${FC} --f77exec=${FC} # --noopt --noarch --debug
  LIB = ${ILOC}
  INC = ${ILOC}
  #F2PY_FLAGS = '-fPIC '

MODULE=pycrtm

#all: wrapped.o\
#	${MODULE}.so

all: ${MODULE}.so

#wrapped.o: wrapped.f90
#	echo ${FORT} whir1
#	echo ${F2PY_COMPILER} whir2
#	${FORT} -c ${FCFLAGS} -L${LIB} -lcrtm  -I${INC} $<
#${MODULE}.so: pycrtm.f90
#	${F2PY} ${F2PY_FLAGS} -c -L${LIB} -lcrtm -lgomp -I${INC} wrapped.o -m ${MODULE}  $< only: wrapForward wrapKmatrix

${MODULE}.so: pycrtm.f90
	${F2PY} ${F2PY_FLAGS} -c -L${LIB} -lcrtm -lgomp -I${INC} -m ${MODULE} $<  
#	mv ${MODULE}*.so ${MODULE}.so 
clean:
	${RM} ${MODULE}*.so
#	${RM} wrapped.o
