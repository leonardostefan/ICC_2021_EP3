    CC     = gcc -g -std=c11
    CFLAGS = 
    LFLAGS = -lm

      PROG = exe3 
      OBJS = odeSolver.o

.PHONY: limpa faxina clean purge all

%.o: %.c %.h 
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

limpa clean:
	@rm -f *~ *.bak

faxina purge:   limpa
	@rm -f *.o core a.out
	@rm -f $(PROG)

