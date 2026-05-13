# Tips to run the programs on a Mac OS

## C++ Programs

- Open Xcode, create a new "Project" select "command inline" give it a name etc.  When all is done a main.cpp is automatically included. Try to run it. If it works then drag a .cpp file from this archive next to the main.cpp. remove the main.cpp. It should run.

- Some of the programs require the file "_kappa.txt". Best is to put it in or next to your project and change the path in the program.

- For Automatic Differentiation you need to put in the same folder as your program file the file "ddouble.hpp".

## Python and Fenics Programs

You may use the grzaphic user interface (GUI) "VS Code" (download it free of charge from the internet) 

## Matlab: use their GUI.

## Fortran

Install the fortran compiler with "homebrew" by typing in the terminal "brew install gfortran"  (Brew is available here https://brew.sh/). Then type 

gfortran myfile.f90 -o myprogram

## FreeFem++

- Download freefem from https://freefem.org/; see the section Assets (for Mac  maybe the version is in the Assets of one version older)

- Follow the instructions. On has to bypass Apple's security rules because freefem is not registered with Apple (it can't be because some shared libraries in Fortran).  Whenever Apple suggests to throw the file away, don't and keep on trying.
  
- Easiest is to use the gui VS Code. Install the "freefem" module (from the module install icon on the left paneel). If your program file ends with .edp or .md you should see a colored syntax and a button to execute the file. If it ends with .md you can also active the right panel and see the math which explains the program 