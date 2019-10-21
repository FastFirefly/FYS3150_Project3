# FYS3150_Project3
Welcome to our FYS3150 Respository for project 3

Here you will find the programs, results and report for Project 3 in FYS3150

## To compile and execute 3a, Gauss Legendre, use: 

Compile:

> g++ -O2 3a.cpp lib.cpp -o 3a.exe -std=c++11

Run: 

> ./3a.exe [Max integration limit] [Amount of mesh points]

Example: 

> g++ -O2 3a.cpp lib.cpp -o 3a.exe -std=c++11
> ./3a.exe 3 20

## To compile and execute 3b, Gauss Laguerre, use: 

Compile:

> g++ -O2 3b.cpp lib.cpp -o 3b.exe -std=c++11

Run: 

> ./3b.exe [Amount of mesh points]

Example: 

> g++ -O2 3b.cpp lib.cpp -o 3b.exe -std=c++11
> ./3b.exe 20


## To compile and execute 3c, Brute Force Monte Carlo, use: 

Compile:

> g++ -O2 3c.cpp -o 3c.exe -std=c++11

Run: 

> ./3c.exe [Max integration limit] [Amount of Integration points]

Example: 

> g++ -O2 3c.cpp -o 3c.exe -std=c++11
> ./3c.exe 3 10000000

## To compile and execute 3d, Importance Sampling Monte Carlo, use: 

Compile:

> g++ -O2 3d.cpp -o 3d.exe -std=c++11 -fopenmp

Run: 

> ./3d.exe [Amount of Integration points]

Example: 

> g++ -O2 3d.cpp -o 3d.exe -std=c++11 -fopenmp
> ./3d.exe 10000000
