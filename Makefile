FC = mpiifort
FLAGS = -c -O3
SUFFIX = 
LF = -warn all
HOME = /home/zxli_1/Lib_90_new

LIBS = $(HOME)/Modules/modules_90.a \
       $(HOME)/MyEis/libeis.a \
       $(HOME)/MyNag/libnag.a \
       $(HOME)/MyLin/liblin.a \
       $(HOME)/Ran/libran.a \
       $(HOME)/Blas/libblas.a

LIBS += -qmkl

all:
	cp $(HOME)/Modules/*.mod . ;\
	make -f Compile FC="$(FC)" LF="$(LF)" FLAGS="$(FLAGS)" LIBS="$(LIBS)" SUFFIX="$(SUFFIX)"

clean: 	
	make -f Compile clean ;\
	rm -f *.mod
