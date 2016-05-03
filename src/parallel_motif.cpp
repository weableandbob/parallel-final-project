/*
Notes:
use combination of MPI_Get_count, MPI_Probe, MPI_recv, and MPI_Isend for message passing
    Probe until a message comes in, get the amount of data, allocate space, recv, operate on received data
May need to use MPI_THREAD_SERIALIZED and have mutexes. This will require main thread to spin on Iprobe instead of blocking on probe
Can't use thread pools for remote threads, otherwise deadlocks may occur
Need to work out pseudocode for parallel DFS
Copy data for each treadDFS call -> less chance of messing up data. If too time/memory expensive, can try reference to single memory block instead
*/

extern "C" {
    #include "C-Thread-Pool/thpool.h"
}
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <list>

#include <mpi.h> //MPI apparently does the extern stuff in its header file

//MPI tags
#define MPIT_RANK_DONE 0
#define MPIT_THREAD_DONE 1
#define MPIT_THREAD_MOVE 2

//Macros
#define CLEANUP_TP \
    thpool_destroy(g_local_threads);

#define CLEAN_EXIT \
    CLEANUP_TP \
    MPI_Finalize(); \
    exit(EXIT_FAILURE);

using namespace std;

/***************************************************************************/
/* Structs *****************************************************************/
/***************************************************************************/
struct graph_node {
    int unique_id; //Unique identifier of the node
    int role; //The role assigned to this node
    set<int> neighbors; //Unique IDs of all nodes with an edge from this node to it
};

struct motif_node {
    int motif_node_index; //Node number in the motif being visited
    int role; //The role the motif node should be
};

struct dfs_data {
    int motif_index; //Index in the global counter array
    struct graph_node; //The graph node being visited on this dfs call
    struct motif_node desired_motif_node; //The motif node id and role that this node should match
    //bool first_invocation; //Whether this is the first time DFS is being run
    list< pair<struct motif_node, struct motif_node> > motif_edges; //List of motif edges left to be evaluated
    set<int> visited_nodes; //Set of all node numbers that have been visited in this path
    map<int, int> motif_to_unique; //Map of motif nodes to unique node IDs
};

/***************************************************************************/
/* Global Variables ********************************************************/
/***************************************************************************/
//MPI globals
int mpi_myrank;
int mpi_commsize;

//Threadpools
threadpool g_local_threads; //Searches that started in this rank
//threadpool g_remote_threads; //Searches that started due to an MPI message

//Graph related variables
int g_ranks_done; //Counter used to keep track of how many nodes are done with the current motif since barriers are not an option
vector< vector<pair<struct motif_node, struct motif_node> > > g_motifs; //A vector of motifs. Each motif is a vector of edges
map<int, struct graph_node> g_local_nodes; //A mapping of unique node identifiers to nodes that are stored in this rank
vector<int> g_vtxdist; //Equivalent to vtxdist from ParMETIS. If vtxdist[i] <= node_id_j < vtxdist[i+1], then the node with unique id j is in rank i
                       //Ordered, so use a binary search

//Synchronization
//pthread_mutex_t* m_counter;
//pthread_mutex_t* m_mpi_lock;
//map<pthread_t, sem_t*> g_locked_threads;

/***************************************************************************/
/* threadDFS ***************************************************************/
/***************************************************************************/
//Depth first search meant to be started in a separate thread
/*Pseudocode:
    if current node's role != desired role: (only want to visit nodes of the correct role)
        return
    if current node has not yet been visited:
        add to visited nodes
        add mapping of motif node id to unique id
    if current node's motif node id does not match the desired id: (takes care of erroneously re-visiting nodes when not in a cycle)
        return
    (From here on, we know that we are at a valid node)
    if no edges left in motif: (this node was the last one needed)
        increment count for motif
        return
    check next edge in motif
    if next edge starts at a previously visited node:
        lookup node
        if node is local:
            call threadDFS on that node
            return
        else:
            transfer data to applicable rank and wait
            return
    (next edge starts at the current node)
    Pop the next edge and store it
    if edge goes to a previously visited node: (Optimization for closing a cycle, as we already know which node we need to go to)
        if that node is a neighbor:
            if neighbor is local:
                call theadDFS on that node and wait
                return
            else:
                transfer data to applicable rank and wait
                return
    (Know that the node we're looking for has not yet been visited)
    for each neighbor:
        if neighbor has been visited:
            continue
        if neighbor is local:
            call threadDFS on that node and wait
            return
        else:
            transfer data to applicable rank and wait
            return*/
void* threadDFS(void* d){
    struct dfs_data* data = (struct dfs_data*)d;
    cout << "DFS run on motif " << data->motif_index << endl;

    return NULL;
}

/***************************************************************************/
/* threadDispatcher *******************************************************/
/***************************************************************************/
//Thread function that walks the graph and searches for the given motif
void* threadDispatcher(void* motif_index){
    //Visit each node in the local graph and search for motifs
    for(int i = 0; i < 10; i++){

        struct dfs_data* data = new struct dfs_data;
        data->motif_index = *(int*)motif_index;
        data->first_invocation = true;
        thpool_add_work(g_local_threads, threadDFS, (void*)data);
    }
    //Wait for all locally dispatched threads to finish
    thpool_wait(g_local_threads);
    //Wait for all ranks to finish then notify own rank
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Send(NULL, 0, MPI_INT, mpi_myrank, MPIT_RANK_DONE, MPI_COMM_WORLD);
    return NULL;
}

/***************************************************************************/
/* main ********************************************************************/
/***************************************************************************/

int main(int argc, char* argv[]){
    int rc;

    //Init MPI
    int mpi_thread_support;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_thread_support);
    if(mpi_thread_support != MPI_THREAD_MULTIPLE){
        cout << "Did not receive the requested thread support" << endl;
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

    //Parse inputs
    if(argc != 4){
        cout << "Usage is " << argv[0] << " graph_file motif_file num_local_threads" << endl;
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    int num_local_threads = atoi(argv[3]);
    //int num_remote_threads = atoi(argv[4]);

    //Initialize thread pools
    g_local_threads = thpool_init(num_local_threads);
    //g_remote_threads = thpool_init(num_remote_threads);

    //Test code: create dummy motifs
    for(unsigned int i = 0; i < 5; i++){
        vector< pair<int, int> > motif;
        for(unsigned int j = 0; j < 5; j++){
            motif.push_back(make_pair(j, j));
        }
        g_motifs.push_back(motif);
    }
    //End test code

    //Search for all motifs
    MPI_Status status;
    int* motif_index = new int;
    for(unsigned int i = 0; i < g_motifs.size(); i++){
        //Create dispatcher thread for the motif
        pthread_t dispatcher;        
        *motif_index = i;
        rc = pthread_create(&dispatcher, NULL, threadDispatcher, motif_index);
        if(rc != 0){
            cout << "Failed to create dispatcher thread" << endl;
            CLEAN_EXIT
        }
        rc = pthread_detach(dispatcher);
        if(rc != 0){
            cout << "Failed to detached dispatcher thread" << endl;
            CLEAN_EXIT
        }

        //Wait for any messages
        while(true){
            rc = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(rc != MPI_SUCCESS){
                cout << "Failed to probe with error code " << rc << endl;
                CLEAN_EXIT
            }
            //Check the tags -> Avoid switch/case since we need to break out of infinite loop
            if(status.MPI_TAG == MPIT_RANK_DONE){
                //Clear message out of buffer
                MPI_Recv(NULL, 0, MPI_INT, mpi_myrank, MPIT_RANK_DONE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                break;
            }
            else if(status.MPI_TAG == MPIT_THREAD_MOVE){
                //Create a new job in remote pool
            }
            else if(status.MPI_TAG == MPIT_THREAD_DONE){
                //Unlock the specified thread
            }
            else{
                cout << "Received unknown tag " << status.MPI_TAG << endl;
                CLEAN_EXIT
            }
        }
        //TODO: gather counts to rank 0
        MPI_Barrier(MPI_COMM_WORLD);
    }
    delete motif_index;
    CLEANUP_TP

    MPI_Barrier(MPI_COMM_WORLD);
    rc = MPI_Finalize();
    cout << "Finalize done " << rc << endl;
    return 0;

    //Begin test code
    /*MPI_Barrier( MPI_COMM_WORLD );
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
    thpool_destroy(tp);*/
    //End test code
}