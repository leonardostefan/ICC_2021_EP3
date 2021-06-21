    CC     = gcc  -std=c11 -g -O3
    CFLAGS = 
    LFLAGS = -lm

      PROG = exe3 
      OBJS = odeSolver.o\
	  		 utils.o

.PHONY: limpa faxina clean purge all

%.o: %.c %.h utils.h Makefile 
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

limpa clean:
	@rm -f *~ *.bak

faxina purge:   limpa
	@rm -f *.o core a.out
	@rm -f $(PROG)

