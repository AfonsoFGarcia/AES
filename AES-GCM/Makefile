rm=/bin/rm -f
CC=cc
DEFS=
INCLUDES=-I.
LIBS=

DEFINES= $(INCLUDES) $(DEFS)
CFLAGS= -std=c99 $(DEFINES) -O2 -fomit-frame-pointer -funroll-loops

all: aes128gcm_driver

aes128gcm_driver: aes128gcm_driver.c aes128e.o aes128gcm.o
	$(CC) $(CFLAGS) -o aes128gcm_driver aes128gcm.o aes128e.o aes128gcm_driver.c 


aes128e.o: ../AES-128/aes128e.c ../AES-128/aes128e.h
	$(CC) $(CFLAGS) -c ../AES-128/aes128e.c $(LIBS)

aes128gcm.o: aes128gcm.c aes128gcm.h
	$(CC) $(CFLAGS) -c aes128gcm.c $(LIBS) 

clean:
	$(rm) aes128e.o aes128gcm_driver *.o core *~

