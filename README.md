* * *

### CICM-TOOLS, version 1.5 for Pure Data

* * *


### About:

 The CICM-Tools are a set of 3 objects allowing you to virtually position a sound source in a multiphonic space. They are:

- [ambicube~] - 3D spatialisation by ambisonic B format
- [ambipan~] - 2D spatialisation by ambisonic B format
- [vbapan~] - 2D spatialisation by Vector Base Amplitude Panning

This project resides at: https://github.com/porres/pd-cicm-tools

Note: CICM-Tools is abandoned and unsupported. The archived version for MAX is found at <https://github.com/CICM/CicmTools>. The Pd download of the original sources and binaries (32 bits only) can be found at <http://cicm.mshparisnord.org/>. 

The version 1.5 for Pure Data was released on july 14th 2004

The port to this new repository was made by Alexandre Torres Porres and aims at "resurecting" this code and also provide binaries for 64 bits as well.

* * *


### Credits :

Copyright (C) 2003-2004 Rémi Mignot, MSH Paris Nord (Maison des Sciences de l'Homme), with Anne Sèdes, Benoît Courribet and Jean-Baptiste Thiebaut, CICM, Paris 8 university, ACI Jeunes Chercheurs "Espaces Sonores".

* * *

### Licence :

The CICM-Tools are under the terms of the GNU LIBRARY GENERAL PUBLIC LICENSE Version 2: http://www.gnu.org/copyleft/gpl.html

* * *

#### Building CICM-Tools for Pd Vanilla:

This project relies on the build system called "pd-lib-builder" by Katja Vetter (see: <https://github.com/pure-data/pd-lib-builder>). PdLibBuilder tries to find the Pd source directory at several common locations, but when this fails, you have to specify the path yourself using the pdincludepath variable. Example:

<pre>make pdincludepath=~/pd-0.51-4/src/  (for Windows/MinGW add 'pdbinpath=~/pd-0.51-4/bin/)</pre>

* Installing with pdlibbuilder

Go to the pd-cicm-tools folder and use "objectsdir" to set a relative path for your build, something like:

<pre>make install objectsdir=../cicm-tools-build</pre>

Then move it to your preferred install folder for Pd and add it to the path.

Cross compiling is also possible with something like this

<pre>make CC=arm-linux-gnueabihf-gcc target.arch=arm7l install objectsdir=../</pre>

* * *

#### 