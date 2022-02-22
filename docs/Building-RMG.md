RMG uses [Cmake](http://www.cmake.org) to manage the build process.[Cmake](http://www.cmake.org) is open source software and is available for a wide range of platforms. Building RMG is nearly automatic on most recent Linux distributions as long as the required packages are installed. Building on various cluster and supercomputer platforms often requires some additional work to get the environment properly configured. We will attempt to provide up to date instructions for various platforms at the the [issue tracker](https://github.com/RMGDFT/rmgdft/issues?q=is%3Aopen+is%3Aissue+label%3A%22Build+instructions%22).

### Supported platforms
* Most modern linux distributions on X86_64 hardware platforms.
* Power PC platforms.
* Cray XK/XE.
* Microsoft windows using the Windows Subsystem for Linux layer.

### Minimum prerequisites
* Cmake v3.1 (earlier versions may work but have not been extensively tested)
* Hardware platform with a 64 bit architecture.
* C++ compiler that supports the C++11 standard. GCC versions 4.8.1 and higher are known to work.
* Boost libraries. Versions 1.5.3 and higher.
* Blas libraries.
* Message passing (MPI) library that supports MPI_THREAD_SERIALIZED.

### Optional packages
* Cuda GPU computing libraries version 9.0 or higher.
* Openbabel chemistry toolbox.
* PLplot scientific plotting package.

### Build instructions
Obtain the code from the RMG git [repository](https://github.com/RMGDFT/rmgdft)
If you have git installed on your local machine use <BR><BR>
`git clone https://github.com/RMGDFT/rmgdft.git`<BR><BR>
this will download the most recent code (not necessarily the most stable code).<BR><BR>
After downloading/unpacking the code or cloning it cd into the top level directory and execute the following instructions from the command line.<BR><BR>
`mkdir build`<BR>
`cd build`<BR>
`cmake ..`<BR><BR>
If you are trying to build a GPU enabled version of RMG change the last line to<BR><BR>
`cmake -DRMG_GPU_ENABLED=1 ..`<BR><BR>
If cmake completes without errors you can now try building RMG using <BR>                             
`make -jN rmg-cpu`<BR><BR>
or if you have used -DRMG_GPU_ENABLED=1<BR><BR>
`make -jN rmg-gpu`<BR><BR>
where N is the number of CPU cores to use in a parallel build.

### Debug builds
To generate a debug build edit the main CMakeList.txt file in the root directory of the distribution and change<BR><BR>
` set(CMAKE_BUILD_TYPE Release)`<BR>
`#set(CMAKE_BUILD_TYPE Debug)`<BR>
to<BR>
`#set(CMAKE_BUILD_TYPE Release)`<BR>
`set(CMAKE_BUILD_TYPE Debug)`<BR>

### Testing
After you have built the rmg-cpu executable you can run system level tests by executing.<BR>
`make test`<BR>
In order to successfully complete all of the tests your machine should have at least 16 GBytes of RAM. The test script will try to use as many CPU cores as your test machine has up to a maximum of 24. On a 24 core Ryzen threadripper it takes around 8 minutes to complete all of the tests.
 
### Troubleshooting
If your platform supports modules you may not have the correct set of modules loaded. Try the command `module avail` to see if your system supports modules and what is available.
If you run into problems you can ask for help at the [issue tracker.](https://github.com/RMGDFT/rmgdft/issues?q=is%3Aopen+is%3Aissue+label%3A%22Build+instructions%22)<BR>
