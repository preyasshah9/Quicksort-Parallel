#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#define ARRAY_SIZE 1048576
#define PROC_NUM 16
#define LOCAL_ARRAY_SIZE 65536
#define PIVOTMETHOD 1
//Pivot Method : 1 -> Using Medians
//0 -> Random Pivoting

//Swapping function
//Used for finding medians and parallel partitioning
void swap(int *a, int *b)
{
    int t;
    t = *a; *a = *b; *b = t;
}

//Compare Function
//Used for quicksort in-built function
int compare (const void *a, const void *b)
{
    return (*(int *) a - *(int *) b);
}

//Partition function
//Used for finding the medians
int partition( int a[], int l, int r) {
    int pivot, i, j, t;
    if(l > r)
    {
        t = rand() % (l - r);
        pivot = a[t];
        swap(&a[t], &a[l]);
    }
    else
    {
        pivot = a[l];
    }
    i = l; j = r+1;
    while( 1)
    {
        do ++i; while( a[i] <= pivot && i <= r );
        do --j; while( a[j] > pivot );
        if( i >= j ) break;
        swap(&a[i], &a[j]);
        //t = a[i]; a[i] = a[j]; a[j] = t;
    }
    swap(&a[l], &a[j]);
    //t = a[l]; a[l] = a[j]; a[j] = t;
    return j;
}

//KthSmallest, function used to find the medians
//Recursive Function
int kSmallest(int A[], int start, int end, int k) {
    int pivotPosition = partition(A, start, end);
    if(pivotPosition == k) {
        return A[k];
    }
    else {
        if(pivotPosition > k) {
            return kSmallest(A, start, pivotPosition - 1, k);
        }
        else {
            return kSmallest(A, pivotPosition + 1, end, k);
        }
    }
}

//Calculate the median of an array
int median(int A[], int numOfElements) {
    //If total number of elements are even, median is avg of n/2th and
    //n/2 - 1th element
    if(numOfElements % 2 == 0) {
        return ((kSmallest(A, 0, numOfElements - 1, numOfElements/2 - 1) +
                 kSmallest(A, 0, numOfElements - 1, numOfElements/2)) / 2);
    }
    else {
        return kSmallest(A, 0, numOfElements - 1, numOfElements/2);
    }
}

int parallel_Partition(int A[], int size, int pivot) {
    int i = 0, j = size;
    while(1) {
        do ++i; while( A[i] <= pivot && i <= size );
        do --j; while( A[j] > pivot );
        if( i >= j ) break;
        swap(&A[i], &A[j]);
    }
    swap(&A[0], &A[j]);
    if(A[j] > pivot) {
        j--;
    }
    return j;
}

void main(int argc, char **argv) {
    int i, j, node;
    MPI_Comm mpiComm;
    int groupLen, arraySize;
    int *totalA;
    int *A, *arrSizeMap, *arrSizePrefixSums, *procArraySize, *procArrayOffset;
    MPI_Init(&argc, &argv);
    clock_t begin, end;
    A = malloc(sizeof(int) * LOCAL_ARRAY_SIZE);
    groupLen = PROC_NUM;
    arraySize = LOCAL_ARRAY_SIZE;
    mpiComm = MPI_COMM_WORLD;
    //Get number of nodes using MPI_Comm_rank Function
    //Get the rank of the current process
    MPI_Comm_rank(mpiComm, &node);

    //Initialize Input Arrays and sends data to different processes
    if(node == 0) {
        totalA = malloc(sizeof(int) * ARRAY_SIZE);
        for (i = 0; i < ARRAY_SIZE; i++) {
            totalA[i] = rand() % 1000000;
        }
    }
    MPI_Scatter(totalA, LOCAL_ARRAY_SIZE, MPI_INT, A, LOCAL_ARRAY_SIZE, MPI_INT, 0, mpiComm);
    //Iterative Quicksort Algorithm
    //Iteration continues until each group of processes have more than one processors
    begin = clock();
    while(groupLen > 1) {
        //Local Variables
        int medianVal;
        int *recvArray;
        int *medianArray,*prefixSumHigherArray, *prefixSumLowerArray;
        int *recvDsplmnt, *sendDsplmnt;
        int *sendCnts, *recvCnts;
        int pivotLocation;
        int numHighElements, numHighProcesses;
        int numLowElements, numLowProcesses;
        int totalElements;
        int nodeCnt, prSumFromNeighbour, arraySizeFromNeighbor;
        int prefixSumLowVal, prefixSumHighVal;

        //Pivoting is done using two methods
        //First, the median is selected from each processor
        //And transmitted to the root process
        if(PIVOTMETHOD == 1) {
            medianVal = median(A, arraySize);
        }
        else {
            //Second method is random sample is send to the root process
            medianVal = A[rand() % arraySize];
        }
        if(node == 0) {
            medianArray = malloc(sizeof(int) * groupLen);
        }
        //Gather all median values into the median array
        MPI_Gather(&medianVal, 1, MPI_INT, medianArray, 1, MPI_INT, 0, mpiComm);
        if(node == 0) {
            //Find the median of pivots sent from each processes into root process
            if(PIVOTMETHOD == 1) {
            	medianVal = median(medianArray, groupLen);
            }
            else {
                medianVal = medianArray[rand() % groupLen];
            }
        }
        //Once a global pivot is selected, broadcast it to the all processes in the group
        MPI_Bcast(&medianVal, 1, MPI_INT, 0, mpiComm);
        //Get the partitioning done at each individual process
        //Based on the global pivot
        //The variable pivotLocation is helpful to further partition the processes
        pivotLocation = parallel_Partition(A, arraySize, medianVal);
        ++pivotLocation; //pivotLocation now contains total number of elements < pivot
        //Root node will perform the prefix-sums and send them to all processes
        prefixSumLowerArray         = malloc(sizeof(int) * (groupLen + 1));
        prefixSumHigherArray        = malloc(sizeof(int) * (groupLen + 1));
        if(node == 0) {
            //Define the prefixSumArray in each send the array to all other processes
            prefixSumLowerArray[0]      = 0;
            prefixSumHigherArray[0]     = 0;
            nodeCnt                     = 1;
            prefixSumLowVal             = pivotLocation;
            prefixSumHighVal            = arraySize - pivotLocation;
            prefixSumLowerArray[nodeCnt]= prefixSumLowVal;
            prefixSumHigherArray[nodeCnt] = prefixSumHighVal;
            while(nodeCnt < groupLen) {
                MPI_Recv(&prSumFromNeighbour, 1, MPI_INT, nodeCnt, 0, mpiComm, MPI_STATUS_IGNORE);
                MPI_Recv(&arraySizeFromNeighbor, 1, MPI_INT, nodeCnt, 0, mpiComm, MPI_STATUS_IGNORE);
                ++nodeCnt;
                prefixSumLowVal += prSumFromNeighbour;
                prefixSumLowerArray[nodeCnt]  = prefixSumLowVal;
                prefixSumHighVal += arraySizeFromNeighbor - prSumFromNeighbour;
                prefixSumHigherArray[nodeCnt] = prefixSumHighVal;
            }
        }
        else {
            MPI_Send(&pivotLocation, 1, MPI_INT, 0, 0, mpiComm);
            MPI_Send(&arraySize, 1, MPI_INT, 0, 0, mpiComm);
        }
        //Once the prefix sum arrays are calculated, broadcast the prefix-sum arrays
        //to All processes in the group
        MPI_Bcast(prefixSumLowerArray, groupLen + 1, MPI_INT, 0, mpiComm);
        MPI_Bcast(prefixSumHigherArray, groupLen + 1, MPI_INT, 0, mpiComm);

        //Based on the pivoting, now determine how the next partition is supposed to be
        //arrSizeMap and arrSizePrefixSums are used to calculate the data that needs to
        //be sent during MPI_alltoallv.
        numHighElements = prefixSumHigherArray[groupLen];
        numLowElements = prefixSumLowerArray[groupLen];
        totalElements = numHighElements + numLowElements;
        arrSizeMap          = malloc(sizeof(int) * groupLen);
        arrSizePrefixSums   = malloc(sizeof(int) * (groupLen + 1));
        arrSizePrefixSums[0] = 0;
        numLowProcesses = (int)round((float)(numLowElements * groupLen)/(float)totalElements);
        if(numLowProcesses == groupLen) {
            numLowProcesses--;
        }
        numHighProcesses = groupLen - numLowProcesses;

        //Populate Array arrSizeMap with the array sizes after partitioning
        //This array will make the calculations of partitioning really easy
        for(i = 0; i < groupLen; i++) {
            if(i < numLowProcesses) {
                if(i == numLowProcesses - 1) {
                    j = i - 1;
                    while(j >= 0) {
                        numLowElements -= arrSizeMap[j];
                        --j;
                    }
                    arrSizeMap[i] = numLowElements;
                }
                else {
                    arrSizeMap[i] = numLowElements/numLowProcesses;
                }
            }
            else {
                if(i == groupLen-1) {
                    j = i - 1;
                    while(j >= numLowProcesses) {
                        numHighElements -= arrSizeMap[j];
                        --j;
                    }
                    arrSizeMap[i] = numHighElements;
                }
                else {
                    arrSizeMap[i] = numHighElements/numHighProcesses;
                }
            }
        }

        //Figure out logic for splitting the processes
        //First, initialize the arrays to be used in MPI_AlltoAllv.
        sendCnts    = malloc(sizeof(int) * groupLen);
        recvCnts    = malloc(sizeof(int) * groupLen);
        sendDsplmnt = malloc(sizeof(int) * groupLen);
        recvDsplmnt = malloc(sizeof(int) * groupLen);
        //Initialize the sendCnt to all 0s.
        //Once arrSizeMap is calculated, calculate the prefixSums
        //Used to populate the sendCnts array
        for(i = 0; i < groupLen; i++) {
            sendCnts[i] = 0;
            arrSizePrefixSums[i + 1] = arrSizeMap[i] + arrSizePrefixSums[i];
        }
        //Now figure out the sendCnt values. populate as per prefix Sums and arrSize map
        //For lower processes
        for(i = 0; i < numLowProcesses; i++) {
            if(prefixSumLowerArray[node] > arrSizePrefixSums[i + 1]) {
                continue;
            }
            else {
                if(prefixSumLowerArray[node + 1] > arrSizePrefixSums[i+1]) {
                    sendCnts[i] = arrSizePrefixSums[i + 1] - prefixSumLowerArray[node];
                    sendCnts[i + 1] = prefixSumLowerArray[node + 1] - arrSizePrefixSums[i + 1];
                    break;
                }
                else {
                    sendCnts[i] = prefixSumLowerArray[node + 1] - prefixSumLowerArray[node];
                    break;
                }
            }
        }
        //I did not handle it very well. But for now, let's just say, prefixSumHigherArray is a bit
        //differently calculated. I just wanted to keep both of these logics separately. They can be
        //Combined
        for(i = numLowProcesses; i < groupLen; i++) {
            if((prefixSumHigherArray[node] + prefixSumLowerArray[groupLen]) > arrSizePrefixSums[i + 1]) {
                continue;
            }
            else {
                if((prefixSumHigherArray[node + 1] + prefixSumLowerArray[groupLen]) > arrSizePrefixSums[i+1]) {
                    sendCnts[i] = arrSizePrefixSums[i + 1] - prefixSumHigherArray[node] - prefixSumLowerArray[groupLen];
                    sendCnts[i + 1] = prefixSumHigherArray[node + 1] + prefixSumLowerArray[groupLen] - arrSizePrefixSums[i + 1];
                    break;
                }
                else {
                    sendCnts[i] = prefixSumHigherArray[node + 1] - prefixSumHigherArray[node];
                    break;
                }
            }
        }
        //Once all sendCnts are populated, now do alltoall communication and populate the recvCnt Array
        MPI_Alltoall(sendCnts, 1, MPI_INT, recvCnts, 1, MPI_INT, mpiComm);

        //Based on the sendCnts and recvCnts, now calculate the send and recv displacements
        sendDsplmnt[0] = 0;
        recvDsplmnt[0] = 0;
        for(i = 0; i < groupLen-1; i++) {
            sendDsplmnt[i + 1] = sendCnts[i] + sendDsplmnt[i];
            recvDsplmnt[i + 1] = recvCnts[i] + recvDsplmnt[i];
        }
        //Calculate the size of the recvBuffer
        int sum = 0;
        for(i = 0; i < groupLen; i++) {
            sum += recvCnts[i]; 
        }
        recvArray = malloc(sizeof(int) * sum);
        MPI_Alltoallv(A, sendCnts, sendDsplmnt, MPI_INT, recvArray,recvCnts, recvDsplmnt, MPI_INT, mpiComm);
        //Once all data is available in the recvBuffer, we don't need A anymore
        free(A);
        //Reallocate the memory for A, and assign elements from recvArray for the next iteration
        A = malloc(sizeof(int) * sum);
        arraySize = sum;
        for(i = 0; i < sum; i++) {
            A[i] = recvArray[i];
        }

        //Cleanup of all pointers
        free(arrSizeMap);
        free(prefixSumLowerArray);
        free(prefixSumHigherArray);
        free(recvDsplmnt);
        free(sendDsplmnt);
        free(sendCnts);
        free(recvCnts);
        if(node == 0) {
            free(medianArray);
        }
        //Now, partition the processes into two groups
        MPI_Comm_split(mpiComm, node < numLowProcesses, node, &mpiComm);
        //Get the upto date process rank and number of processes in group
        MPI_Comm_size(mpiComm, &groupLen);
        MPI_Comm_rank(mpiComm, &node);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    qsort(A, arraySize, sizeof(int), compare);

    //All of these code is to gather data in process 0 so that it can be printed nicely
    if(node == 0) {
        end = clock();
        printf("Execution Time: %lf\n", (double)(end - begin)/ (double)CLOCKS_PER_SEC);
        procArraySize   = malloc(sizeof(int) * PROC_NUM);
    }
    MPI_Gather(&arraySize, 1, MPI_INT, procArraySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(node == 0) {
        procArrayOffset = malloc(sizeof(int) * PROC_NUM);
        procArrayOffset[0] = 0;
        for(i = 1; i < PROC_NUM; i++) {
            procArrayOffset[i] = procArrayOffset[i-1] + procArraySize[i-1];
        }
    }
    MPI_Gatherv(A, arraySize, MPI_INT, totalA, procArraySize, procArrayOffset, MPI_INT, 0, MPI_COMM_WORLD);
    if(node == 0) {
        for(i = 0; i < ARRAY_SIZE; i++) {
            printf("%d\t", totalA[i]);
        }
    }
    
    free(A);
    if(node == 0) {
        free(totalA);
        free(procArrayOffset);
        free(procArraySize);
    }
    MPI_Finalize();
}
