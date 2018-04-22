# Parralel-MCF
A university project to create a parallelization recipe using the MCF benchmark

## Building
Windows:  

    gcc -g .\*.c -o mcf.exe -fopenmp
## Running
### Running on train input:  
    mcf data/train/input/inp.in  
### Running on reference input:  
    mcf data/ref/input/inp.in
