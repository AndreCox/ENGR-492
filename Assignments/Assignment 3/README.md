# Assignment 3

This project was built with modern C++ on Fedora Linux. However it should compile for other OS’s. I have tested on Fedora 42 and Ubuntu 22.04.3 LTS
## Dependencies 
- CMake
- Make
- gcc & g++
- git

## Linux build dependencies
These can be installed using apt or dnf.

**Fedora/RHEL/CentOS**

```libX11-devel libXrandr-devel libXcursor-devel libXi-devel  mesa-libGL-devel libXxf86vm-devel libXinerama-devel libXfixes-devel systemd-devel freetype-devel libvorbis-devel libogg-devel flac-devel```

**Ubuntu/Debian**

```libx11-dev libxrandr-dev libxcursor-dev libxi-dev libgl1-mesa-dev libxinerama-dev libxfixes-dev libudev-dev libfreetype6-dev libvorbis-dev libogg-dev libflac-dev```


Start by cloning the git repo

```sh
git clone https://github.com/AndreCox/ENGR-492.git
```

Enter the path for the project

```sh
cd ENGR-492/Assignments/Assignment\ 2/
```


Make and enter the build directory

```sh
mkdir build && cd build
```

Generate build files

```sh
cmake ..
```

Once this completes you can compile the code

```sh
make
```

After compilation finishes you can run the code with

```sh
./main
```

