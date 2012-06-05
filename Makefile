#DEBUGFLAG=-debug
OPTIMIZE=-O -inline

CMD = rebuild $< $(OPTIMIZE) $(DEBUGFLAG) -Isrc -Idev -oqobj -ofbin/$@

all: converting analysing generating testing

converting: gmapconvert gmappack
analysing:  gmaphardyweinberg gmapassoc gmapfreq
generating: gmapgenerate gmaprandpheno 
testing:    gmapmulassoc

gmapconvert: src/gmapconvert.d
	$(CMD)

gmapgenerate: src/gmapgenerate.d
	$(CMD)

gmaphardyweinberg: src/gmaphardyweinberg.d
	$(CMD)

gmapassoc: src/gmapassoc.d
	$(CMD)

gmapfreq: src/gmapfreq.d
	$(CMD)

gmappack: src/gmappack.d
	$(CMD)

gmaprandpheno: src/gmaprandpheno.d
	$(CMD)

gmapmulassoc: src/gmapmulassoc.d
	$(CMD)

cleanall: clean
clean:
	rm -f obj/*
