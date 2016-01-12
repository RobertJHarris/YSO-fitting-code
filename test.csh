#!/bin/csh
gcc -I./ -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -fno-strict-aliasing -fno-common -dynamic  -arch x86_64 -g -Os -pipe -fno-common -fno-strict-aliasing -fwrapv -DENABLE_DTRACE -DMACOSX -DNDEBUG -Wall -Wstrict-prototypes -Wshorten-64-to-32 -DNDEBUG -g -fwrapv -Os -Wall -Wstrict-prototypes -DENABLE_DTRACE -c fitsio.c 

gcc -arch x86_64 -c make_vis.c -o make_vis.o

gcc -g -I./ -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -fno-strict-aliasing -fno-common -dynamic -arch x86_64 -g -Os -pipe -fno-common -fno-strict-aliasing -fwrapv -DENABLE_DTRACE -DMACOSX -DNDEBUG -Wall -Wstrict-prototypes -Wshorten-64-to-32 -DNDEBUG -g -fwrapv -Os -Wall -Wstrict-prototypes -DENABLE_DTRACE -L/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config -L/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/ -lpython2.7 -ldl  -lfftw3 -lm make_vis.o fitsio.o create_model.c  -o create_model


./create_model 1024 test.fits test512.fits 11 > new

ds9 after.fits padded_real.fits
