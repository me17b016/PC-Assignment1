# PC-Assignment1
`Find shortest path using Matrix-Matrix Multiplication`
# Plagiarism Self Declaration:
`I, Raj Garg , certify that this assignment is my own work,
and entirely based on my personal study and understanding. I have not copied in part or
whole or otherwise plagiarised the work of other students and/or persons.`
# Configuration
#### Conditions
Let say **P** is the number of processes, then `P must be perfect square`

Let say **N** is the number of nodes , then `(N * N) % P == 0` 

#### Setup
Uncomment the line no. **19, 22, 121**, If the files [Solve.cpp, HelloMPI.cpp] are in the same directory Or mention the `path` in line no. **18, 21, 120**

Compile `HelloMPI.cpp`, and Executable file must be same directory. Or mention the Path as stated above.

Run only `Solve.cpp`, no need to run  `HelloMPI.exe`

`Solve.cpp` will call `HelloMPI.exe` itself by using this command `mpiexec -np [processes] HelloMPI.exe A.txt G.txt`, You may need to change `mpiexec` according to your system in `Solve.cpp` file, line no. **84**.

No need to add `A.txt` and `G.txt`, those file will generate by `Solve.cpp`

#### Output
You should get output like this:

` Enter number of Nodes :`

`16`

`Enter Source vertex:`

`1`

`Enter Destination vertex:`

`8`

`Enter number of process you want to run :16`

`Elapsed time by Rank 0 : 0.0100341000`

`Path length is 2`
