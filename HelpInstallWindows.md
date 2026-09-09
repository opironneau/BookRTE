# Tips to run the programs on a Windows PC

## C++ Programs

- Install MSYS2 from https//www.MSYS2.org. You need to type in the terminal of MSYS2  
 
 pacman -S mingw-w64-ucrt-x86_64-gcc

The Terminal command is g++ myfile.cpp.  If you want an IDE (integrated Development Environment) with a debugger you can try Visual Studio Code but I found its installation requires some help: see https://www.youtube.com/watch?v=5I3WKjDHox0) 

- Some of the programs require the file "_kappa.txt". Best is to put it in or next to your project and change the path in the program.

- For Automatic Differentiation you need to put in the same folder as your program file the file "ddouble.hpp".

## Python and Fenics Programs

You may use VS Code (short for Visual Studio Code);  download it free of charge from   https://code.visualstudio.com/download

## Matlab: 
use their IDE.

## Fortran

Install the fortran compiler with the command in the MSYS2 window

pacman -S mingw-w64-x86_64-gcc-fortran

To compile do

gfortran myfile.f90 -o myprogram

If it doesn't work look into  C:\mysys64\ for another terminal (shell) of MSYS2 called MSYS2 MINGW64 SHELL  with a blue icon.

## FreeFem++

- Download freefem from https://freefem.org/ . One has to bypass Windows' security rules because freefem is not registered with Microsoft (it can't be because some shared libraries are in Fortran).  If the install file is not called 
  
  FreeFEM-4.16-amd64-win64.exe"
  
  it is because Windows changed its name for security.  You must change the name back to the above or at least change its extension to .exe.
  
- Easiest is to use the gui VS Code. Install the "freefem" extension module (from the extension icon on the left panel). If your program file ends with .edp or .md you should see a colored syntax and a button to execute the file. If it ends with .md you can also activate the right panel and see the math which explains the program 