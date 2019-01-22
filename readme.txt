QuickSort Algorithm: main.c
Jobscript and makefiles are attached. 
Put all of these files in your home directory and execute command: 
make.
QuickSort Algorithm as described in book 409-410.
Currently main.c and jobscript is configured as sorting 1 million integrers.
Number of processes 16. If number of processes is changed, then update
the jobscript needs to be updated accordingly.
In code, please update the #define statements as: 
PROC_NUM -> number of processes.
LOCAL_ARRAY_SIZE -> ARRAY_SIZE/PROC_NUM

To switch between pivot methods, change #define statement:
#define PIVOTMETHOD 0/1
0 -> Random Pivoting
1 -> Median Pivoting


