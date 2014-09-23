%module cgp
 %{
 /* Includes the header in the wrapper code */
 #include "../src/cgp.h"
 %}
 
 /* Parse the header file to generate wrappers */
 %include "../src/cgp.h"
