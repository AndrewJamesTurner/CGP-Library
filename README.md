CGP Library
======

A platform independent Cartesian Genetic Programming Library written in C.

Author: Andrew James Turner  
Email: andrew.turner@york.ac.uk  
License: Lesser General Public License (LGPL)  

##How To Use

###Simplest

The simplest method is to add cgp.c and cgp.h to your build path and compile them along with your own files.

For instance #include cgp.h and then compile your program similar too:

    $ gcc yourFiles.c cgp.c cgp.h

###Standard

Compile source into a shared library and install system wide.

####Ubuntu

To generate the shared object (.so) open a terminal in the CGP-Library directory and run:

    $ make so
    
Then copy the generated libcgp.so to /usr/lib by running:

    $ sudo cp libcgp.so /usr/lib
    
Give libcgp.so the necessary permissions by running:

    $ sudo chmod 0755 /usr/lib/libcgp.so

//Then copy cgp.h to /usr/include by running:

//    $ sudo cp ./src/cgp.h /usr/include
    
Finally update the system so it knows about libcgp.so

    $ sudo ldconfig
    
Once CGP-Library has been installed it can be accessed by including cgp.h in your project and linking using the -lcgp flag.

####Windows

####Mac
