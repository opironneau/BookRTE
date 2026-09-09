# Tips to run the programs on Ubunttu-linux

## C++ Programs

- Install VS Code (Visual Studio Code) https://code.visualstudio.com/download and the Microsoft C++ extension. Make sure the gnu compiler and the debugger anre installed. If not
  
  gcc -v
  
  if not installed

sudo apt-get update

sudo apt-get install build-essential gdb


- Some of the programs require the file "_kappa.txt". Best is to put it in or next to your project and change the path in the program.

- For Automatic Differentiation you need to put in the same folder as your program file the file "ddouble.hpp".

## Python and Fenics

Normally python3 is install with the Ubuntu distribution.  So easiest is to install VS Code, the microsoft python extension module.  If it doesn't work see the documentation of VS Code https://code.visualstudio.com/docs/languages/python


## Matlab

Use their IDE.

## Fortran

Install the fortran compiler 

sudo apt-get install gfortran

 Then type 

gfortran myfile.f90 -o myprogram

## FreeFem++

- Download freefem from 
  
  https://freefem.org/
  
  See the section Assets and download the version compatible with your processor.
  
- To run a freefem script, easiest is to use VS Code (see above). Install the "freefem" module (from the extension install icon on the left panel). If your program file ends with .edp or .md you should see a colored syntax and a button to execute the file. If it ends with .md you can also active the right panel and see the math which explains the program .  If you don't like VS Code you can also run it via the Terminal  by typing 
  
  FreeFem++ myfile.edp