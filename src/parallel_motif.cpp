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
#define MPIT_REQUEST_NEIGHBORS 3
#define MPIT_RESPONSE_NEIGHBORS 4

//Macros
#define CLEANUP_TP \
    thpool_destroy(g_local_threads);

#define CLEAN_EXIT \
    CLEANUP_TP \
    MPI_Finalize(); \
    exit(EXIT_FAILURE);

#define CHECK_MUTEX_INIT(value) \
    if(value != 0){ \
        cerr << "Failed to initialize mutex with error code " << value << endl; \
        CLEAN_EXIT \
    }

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
    struct graph_node cur_node; //The graph node being visited on this dfs call
    struct motif_node desired_motif_node; //The motif node id and role that this node should match
    bool first_invocation; //Whether this is the first time DFS is being run
    vector< pair<struct motif_node, struct motif_node> > motif_edges; //List of motif edges left to be evaluated
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

//Graph related variables
int g_ranks_done; //Counter used to keep track of how many nodes are done with the current motif since barriers are not an option
vector< vector<pair<struct motif_node, struct motif_node> > > g_motifs; //A vector of motifs. Each motif is a vector of edges
map<int, struct graph_node> g_local_nodes; //A mapping of unique node identifiers to nodes that are stored in this rank
vector<int> g_vtxdist; //Equivalent to vtxdist from ParMETIS. If vtxdist[i] <= node_id_j < vtxdist[i+1], then the node with unique id j is in rank i
                       //Ordered, so use a binary search

//Synchronization
pthread_mutex_t* m_counter_lock;
pthread_mutex_t* m_mpi_lock;
//map<pthread_t, sem_t*> g_locked_threads;

/***************************************************************************/
/* isValidMotif ************************************************************/
/***************************************************************************/
//Checks a possible motif for excess edges, aka edges between nodes in the motif that are not defined in the motif
//If any such edges exist, the given set of vertices is not a valid motif
/*Pseudocode:
    for each visited node:
        set motif_neighbors
        if node is local:
            copy neighbors that are also in visited set into motif_neighbors
        else:
            send request to applicable rank for node's neighbors that are also in the visited set
            wait for response from main thread
            copy response into motif_neighbors
        for each neighbor in motif_neighbor:
            if the motif does not have an edge from this node to neighbor:
                return false
        return true
*/
bool isValidMotif(dfs_data* data){
    return true;
}

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
        if not valid motif:
            return
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

    //Test code
    cout << "Rank " << mpi_myrank << " DFS run on motif " << data->motif_index << " graph node role " << data->cur_node.unique_id;
    cout << " desired motif node role " << data->desired_motif_node.motif_node_index << " motif edges ";
    for(unsigned int i = 0; i < data->motif_edges.size(); i++){
        cout << "([" << data->motif_edges[i].first.motif_node_index << ", ";
        cout << data->motif_edges[i].first.role << "], [" << data->motif_edges[i].second.motif_node_index;
        cout << ", " << data->motif_edges[i].second.role << "]) ";
    }
    cout << endl;
    //End test code

    //Current node and desired motif node role mismatch
    if(data->cur_node.role != data->desired_motif_node.role){
        return NULL;
    }

    return NULL;
}

/***************************************************************************/
/* threadDispatcher *******************************************************/
/***************************************************************************/
//Thread function that walks the graph and searches for the given motif
void* threadDispatcher(void* motif_index){
    //Visit each node in the local graph and search for motifs
    for(map<int, graph_node>::iterator it = g_local_nodes.begin(); it != g_local_nodes.end(); it++){
        struct dfs_data* data = new struct dfs_data;

        data->motif_index = *(int*)motif_index;
        data->cur_node = it->second;
        data->first_invocation = true;
        data->motif_edges = g_motifs[data->motif_index];
        data->desired_motif_node = data->motif_edges[0].first;

        thpool_add_work(g_local_threads, threadDFS, (void*)data);
        //Test code
        sleep((mpi_myrank+1) * 3);
        //End test code
    }
    //Wait for all locally dispatched threads to finish
    thpool_wait(g_local_threads);

    //Notify all ranks of this rank's completion
    //This would be so much better if we could just MPI_Barrier
    MPI_Request req = MPI_REQUEST_NULL;
    pthread_mutex_lock(m_mpi_lock);
    for(int i = 0; i < mpi_commsize; i++){
        MPI_Isend(NULL, 0, MPI_INT, i, MPIT_RANK_DONE, MPI_COMM_WORLD, &req);
    }    
    pthread_mutex_unlock(m_mpi_lock);
    return NULL;
}

/***************************************************************************/
/* main ********************************************************************/
/***************************************************************************/

int main(int argc, char* argv[]){
    int rc;

    //Init MPI
    int mpi_thread_support;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_support);
    if(mpi_thread_support != MPI_THREAD_SERIALIZED){
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

    //Initialize thread pools
    g_local_threads = thpool_init(num_local_threads);

    //Initialize synchronization objects
    m_counter_lock = new pthread_mutex_t;
    rc = pthread_mutex_init(m_counter_lock, NULL);
    CHECK_MUTEX_INIT(rc)
    m_mpi_lock = new pthread_mutex_t;
    rc = pthread_mutex_init(m_mpi_lock, NULL);
    CHECK_MUTEX_INIT(rc)

    //Test code: create dummy graph
    for(int i = 0; i < 2; i++){
        struct graph_node v;
        v.unique_id = i;
        v.role = 0;
        for(int j = 0; j < 2; j++){
            if(i == j){
                continue;
            }
            v.neighbors.insert(j);
        }
        g_local_nodes[i] = v;
    }
    //End test code

    //Test code: create dummy motifs
    for(unsigned int i = 0; i < 2; i++){
        vector< pair<struct motif_node, struct motif_node> > motif;
        for(unsigned int j = 0; j < 2; j++){
            struct motif_node a;
            a.motif_node_index = i;
            a.role = 0;
            struct motif_node b;
            b.motif_node_index = j;
            b.role = 0;
            motif.push_back(make_pair(a, b));
        }
        g_motifs.push_back(motif);
    }
    //End test code

    //Search for all motifs
    MPI_Status status;
    int flag;
    int* motif_index = new int;
    for(unsigned int i = 0; i < g_motifs.size(); i++){
        g_ranks_done = 0;

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
            //Look for any incoming messages
            pthread_mutex_lock(m_mpi_lock);
            rc = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            pthread_mutex_unlock(m_mpi_lock);
            if(rc != MPI_SUCCESS){
                cout << "Failed to probe with error code " << rc << endl;
                CLEAN_EXIT
            }
            if(!flag){ //No messages in any queues
                continue;
            }

            //Check the tags -> Avoid switch/case since we need to break out of infinite loop
            if(status.MPI_TAG == MPIT_RANK_DONE){
                //Clear message out of buffer and increment number of ranks done
                pthread_mutex_lock(m_mpi_lock);
                MPI_Recv(NULL, 0, MPI_INT, status.MPI_SOURCE, MPIT_RANK_DONE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                pthread_mutex_unlock(m_mpi_lock);
                g_ranks_done++;
                if(g_ranks_done == mpi_commsize){
                    break; //All ranks are done, so break out of infinite loop
                }
            }
            else if(status.MPI_TAG == MPIT_THREAD_MOVE){
                //Spawn and detach a remote thread
            }
            else if(status.MPI_TAG == MPIT_THREAD_DONE){
                //Unlock the specified thread
            }
            else if(status.MPI_TAG == MPIT_REQUEST_NEIGHBORS){
                //Get requested neighbors and reply
            }
            else if(status.MPI_TAG == MPIT_RESPONSE_NEIGHBORS){
                //Copy response into shared buffer and notify relevant thread
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