CFLAGS = -Wall -W -O2
LDFLAGS = -s

CC = gcc
INSTALL = install
INSTALL_DATA = $(INSTALL) -m 444
INSTALL_PROGRAM = $(INSTALL) -m 555

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
mandir = $(prefix)/share/man
man1dir = $(mandir)/man1
manext = .1

BIN = base91

.PHONY: all install check clean

all: $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

base91: cli.o base91.o
	$(CC) $(LDFLAGS) -o $@ $^

install: all
	mkdir -p $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) base91 $(DESTDIR)$(bindir)/base91
	ln -sf base91 $(DESTDIR)$(bindir)/b91dec
	ln -sf base91 $(DESTDIR)$(bindir)/b91enc
	mkdir -p $(DESTDIR)$(man1dir)
	$(INSTALL_DATA) base91.1 $(DESTDIR)$(man1dir)/base91$(manext)
	ln -sf base91$(manext) $(DESTDIR)$(man1dir)/b91dec$(manext)
	ln -sf base91$(manext) $(DESTDIR)$(man1dir)/b91enc$(manext)

check: all
	cd test && $(MAKE)

clean:
	-rm -f *.o $(BIN) core
	cd test && $(MAKE) clean
