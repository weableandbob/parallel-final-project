/*
Notes:
use combination of MPI_Get_count, MPI_Probe, MPI_recv, and MPI_Isend for message passing
    Probe until a message comes in, get the amount of data, allocate space, recv, operate on received data
*/

extern "C" {
    #include "C-Thread-Pool/thpool.h"
}
#include <iostream>
#include <cstdlib>
#include <unistd.h>

#include <mpi.h> //MPI apparently does the extern stuff in its header file

#define MPIT_TEST 0

using namespace std;

/***************************************************************************/
/* Global Variables ********************************************************/
/***************************************************************************/
//MPI globals
int mpi_myrank;
int mpi_commsize;

/***************************************************************************/
/* test_thread *************************************************************/
/***************************************************************************/
//Used for testing the threadpool library
void* test_thread(void* id){
    sleep(*(int*)id);
    cout << "Hello from thread " << *(int*)id << " with actual id " << pthread_self() << endl;
    delete (int*)id;
    return NULL;
}

/***************************************************************************/
/* main ********************************************************************/
/***************************************************************************/

int main(int argc, char* argv[]){
    //Init MPI
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

    //Begin test code
    MPI_Barrier( MPI_COMM_WORLD );
    int message;
    MPI_Request reqs[2];
    reqs[0] = MPI_REQUEST_NULL;
    reqs[1] = MPI_REQUEST_NULL;
    MPI_Status stats[2];
    MPI_Irecv(&message, 1, MPI_INT, (mpi_myrank + 1)%mpi_commsize, MPIT_TEST, MPI_COMM_WORLD, &reqs[0]);
    MPI_Isend(&mpi_myrank, 1, MPI_INT, (mpi_myrank + mpi_commsize - 1)%mpi_commsize, MPIT_TEST, MPI_COMM_WORLD, &reqs[1]);
    MPI_Waitall(2, reqs, stats);
    cout << "Rank " << mpi_myrank << " received message " << message << endl;


    cout << "Test" << endl;
    threadpool tp = thpool_init(4);
    for(int i = 0; i < 12; i++){
        int* id = new int;
        *id = i;
        thpool_add_work(tp, test_thread, (void*)id);
        
    }
    thpool_wait(tp);
    thpool_destroy(tp);
    //End test code

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}