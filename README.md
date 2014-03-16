CGP Library
======

A platform independent Cartesian Genetic Programming Library written in C.

Author: Andrew James Turner  
Email: andrew.turner@york.ac.uk  
License: Lesser General Public License (LGPL)  

##How To Use

###Simplest

The simplest method is to add libCGP.c and libCGP.h to your build path and compile them along side your own files.

###Standard

Compile source into a shared library and install system wide.

####Debian

To generate the shared object (.so) open a terminal in the CGP-Library directory and run:

    $ make so
    
Then copy the generated libCGP.so to /usr/lib by running:

    $ sudo cp libGP.so /usr/lib
    
Give libCGP.so the necessary permissions by running:

    $ sudo chmod 0755 /usr/lib/libGCP.so

Then copy libCGP.h to /usr/include by running:

    $ sudo cp /src/libCGP.h /user/include
    
Finally update they system so it knows about libCGP.so

    $ sudo ldconfig

####Windows

####Mac
