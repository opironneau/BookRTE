# Tips to run the programs on a Mac OS

## C++ Programs

- Open Xcode, create a new "Project" select "command inline" give it a name etc.  When all is done a main.cpp is automatically included. Try to run it. If it works then drag a .cpp file from this archive next to the main.cpp. remove the main.cpp. It should run.

- Some of the programs require the file "_kappa.txt". Best is to put it in or next to your project and change the path in the program.

- For Automatic Differentiation you need to put in the same folder as your program file the file "ddouble.hpp".

## Python and Fenics

You may use the integrated development environment (IDE) "VS Code". Check your cpu (X86 or Arm64) and download it from

https://code.visualstudio.com/download 

## Matlab

Use their IDE.

## Fortran

Install the fortran compiler with "homebrew" by typing in the terminal "brew install gfortran"  (Brew is available here https://brew.sh/). Then type 

gfortran myfile.f90 -o myprogram

## FreeFem++

- Download freefem from 
  
  https://freefem.org/
  
  See the section Assets (for Mac  it may be that the version is in the Assets of one  older release section).

- Follow the instructions. On has to bypass Apple's security rules because freefem is not registered with Apple (it can't be because some shared libraries in Fortran).  Whenever Apple suggests to throw the file away, don't! see the "security" section in the system parameters (below the apple icon in the very top left menu bar) and keep on trying.
  
- To run a freefem script, easiest is to use VS Code (see above). Install the "freefem" module (from the extension install icon on the left panel). If your program file ends with .edp or .md you should see a colored syntax and a button to execute the file. If it ends with .md you can also active the right panel and see the math which explains the program .  If you don't like VS Code you can edit your .edp file with any editor and drag it to the FreeFem icon (in the application folder) to run it.  You can also execute it via the Terminal (in Applications>Utilities folder) by typing 
  
  FreeFem++ myfile.edp