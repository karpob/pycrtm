module load other/comp/gcc-7.2
source ~/bin/go3.csh 
setenv FORT gfortran
setenv FCFLAGS '-fimplicit-none -ffree-form -fopenmp -fno-second-underscore -fPIC -frecord-marker=4 -std=f2008'
setenv ILOC /discover/nobackup/bkarpowi/rt/crtm_2.3.0/
setenv F2PY_COMPILER gnu95
