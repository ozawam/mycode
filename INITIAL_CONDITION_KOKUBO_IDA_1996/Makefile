export OPENGL=0
#export MPI=1
export CC=gcc
INCDIR  = src_original/
EXEFILE = REBOUND

include $(INCDIR)Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I$(INCDIR) -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lrebound $(LIB) -o $(EXEFILE)
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(INCDIR)
	@-rm -f librebound.so
	@ln -s $(INCDIR)librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(INCDIR) clean
	@echo "Cleaning up local directory ..."
	@-rm -vf $(EXEFILE) 
