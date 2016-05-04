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

#define END_DFS(d) \
    delete d; \
    return NULL;

/*#define COPY_STRUCT_FOR_LOCAL(orig, copy) \
    copy->motif_index = orig->motif_index; \
    copy->cur_node = g_local_nodes[next_node]; \
    copy->desired_motif_node.motif_node_index = it->first; \
    copy->desired_motif_node.role = copy->cur_node.role; \
    copy->first_invocation = false; \
    copy->motif_edges = orig->motif_edges; \
    copy->visited_nodes = orig->visited_nodes; \
    copy->motif_to_unique = orig->motif_to_unique;*/

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

//Graph related variables
int g_ranks_done; //Counter used to keep track of how many nodes are done with the current motif since barriers are not an option
vector< list<pair<struct motif_node, struct motif_node> > > g_motifs; //A vector of motifs. Each motif is a list of edges
map<int, struct graph_node> g_local_nodes; //A mapping of unique node identifiers to nodes that are stored in this rank
vector<int> g_vtxdist; //Equivalent to vtxdist from ParMETIS. If vtxdist[i] <= node_id_j < vtxdist[i+1], then the node with unique id j is in rank i
                       //Ordered, so use a binary search
vector<int> g_motif_counts; //Counts of each motif that were found to end in this local graph

//Synchronization
pthread_mutex_t* m_counter_lock;
pthread_mutex_t* m_mpi_lock;
//map<pthread_t, sem_t*> g_locked_threads;

/***************************************************************************/
/* printGraph **************************************************************/
/***************************************************************************/
void printGraph(){
    for(map<int, struct graph_node>::iterator iter = g_local_nodes.begin();
            iter != g_local_nodes.end(); iter++){
        cout << "Node id " << iter->second.unique_id << " role " << iter->second.role << endl;
        cout << "Neighbors: ";
        for(set<int>::iterator n = iter->second.neighbors.begin(); 
                n != iter->second.neighbors.end(); n++){
            cout << *n << " ";
        }
        cout << endl;
    }
}

/***************************************************************************/
/* getRankForNode **********************************************************/
/***************************************************************************/
//Returns the MPI rank storing the node with the given unique id using a binary search
int getRankForNode(int id){
    if(g_vtxdist[g_vtxdist.size() - 1] <= id){
        cerr << "Requested to find a node not in the graph" << endl;
        CLEAN_EXIT
    }

    int left = 0;
    int right = g_vtxdist.size() - 2;
    int m;
    while(true){
        //Don't need to check if left > right since the format of vtxdist guarantees a node will be inside
        m = (left + right) / 2;
        if(m == (int)(g_vtxdist.size() - 1)){
            m -= 1;
        }
        if(g_vtxdist[m] <= id && g_vtxdist[m+1] > id){
            return m;
        }
        else if(g_vtxdist[m] <= id){
            left = m + 1;
        }
        else{
            right = m - 1;
        }
    }
}

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
    (void)data;
    return true;
}

/***************************************************************************/
/* threadDFS ***************************************************************/
/***************************************************************************/
//Depth first search meant to be started in a separate thread
//The function is responsible for deleting any data it is given
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
    cout << "Rank " << mpi_myrank << " DFS run on motif " << data->motif_index << " graph node " << data->cur_node.unique_id;
    cout << " desired motif node role " << data->desired_motif_node.motif_node_index << " motif edges ";
    for(list<pair<struct motif_node, struct motif_node> >::iterator i = data->motif_edges.begin();
            i != data->motif_edges.end(); i++){
        cout << "([" << i->first.motif_node_index << ", ";
        cout << i->first.role << "], [" << i->second.motif_node_index;
        cout << ", " << i->second.role << "]) ";
    }
    cout << endl;
    //End test code

    //Ensure we're at a valid node for the motif
    //Current node and desired motif node role mismatch
    if(data->cur_node.role != data->desired_motif_node.role){
        cout << "mismatched role" << endl;
        END_DFS(data)
    }
    //Current node is unvisited
    if(data->visited_nodes.find(data->cur_node.unique_id) == data->visited_nodes.end()){
        data->visited_nodes.insert(data->cur_node.unique_id);
        data->motif_to_unique[data->desired_motif_node.motif_node_index] = data->cur_node.unique_id;
    }
    //Check whether the node's motif id is the same as the desired id
    if(data->motif_to_unique[data->desired_motif_node.motif_node_index] != data->cur_node.unique_id){
        cout << "mismatch id" << endl;
        END_DFS(data)
    }

    cout << "valid node" << endl;
    //Now know that the node we're at is valid
    //Check if this was the last node needed to complete motif
    if(data->motif_edges.empty()){
        if(isValidMotif(data)){
            pthread_mutex_lock(m_counter_lock);
            g_motif_counts[data->motif_index] = g_motif_counts[data->motif_index] + 1;
            pthread_mutex_unlock(m_counter_lock);
        }
        END_DFS(data);
    }

    //Check whether the next edge starts at this node or not
    int source_motif_node = data->motif_edges.front().first.motif_node_index;
    map<int, int>::iterator it = data->motif_to_unique.find(source_motif_node);
    //Does not support motifs that have more than node that cannot be reached by any other node
    //E.g. 0 -> 1 <- 2 is invalid since 0 and 2 cannot be reached unless starting the search at them
    if(it == data->motif_to_unique.end()){
        cerr << "Was given a motif that is not supported by this motif representation format" << endl;
        CLEAN_EXIT
    }

    cout << "valid motif representation" << endl;
    //next_node is a unique node identifier
    int source_unique_node = it->second;
    int rank_with_node;
    //Need to hop to previously visited node before continuing
    if(source_unique_node != data->cur_node.unique_id){
        //Determine if the visited node is in this rank or not
        rank_with_node = getRankForNode(source_unique_node);
        //Node is local, so simply jump
        if(rank_with_node == mpi_myrank){
            cout << "jumping" << endl;
            struct dfs_data* input = new dfs_data;
            //Copy data into new struct to pass to new invocation of threadDFS
            input->motif_index = data->motif_index;
            input->cur_node = g_local_nodes[source_unique_node];
            //Jumping, so the motif_node_index should be the same as the previously visited node
            input->desired_motif_node.motif_node_index = it->first;
            //Role should be the same as the previously visited node, as well
            input->desired_motif_node.role = input->cur_node.role;
            input->first_invocation = false;
            input->motif_edges = data->motif_edges;
            input->visited_nodes = data->visited_nodes;
            input->motif_to_unique = data->motif_to_unique;

            threadDFS((void*)input);
            END_DFS(data)
        }
        //Node is in another rank, so send data to that rank and continue
        else{

        }
    }

    cout << "continues from current node" << endl;
    //Otherwise continues from this node
    pair<struct motif_node, struct motif_node> next_edge = data->motif_edges.front();
    //Pop the edge so we don't keep calling dfs on the same edge
    data->motif_edges.pop_front();

    //Determine if we have visited the next node before
    int dest_motif_node = next_edge.second.motif_node_index;
    it = data->motif_to_unique.find(dest_motif_node);

    //The next node was previously visited
    if(it != data->motif_to_unique.end()){
        int dest_unique_node = it->second;
        //If the node is not a neighbor, know we can immediately return
        if(data->cur_node.neighbors.find(dest_unique_node) == data->cur_node.neighbors.end()){
            cout << "previously visited node not neighbor" << endl;
            END_DFS(data)
        }
        //Determine if the visited node is in this rank or not
        rank_with_node = getRankForNode(dest_unique_node);
        //Node is local, so simply call from same thread
        if(rank_with_node == mpi_myrank){
            cout << "previously visited" << endl;
            struct dfs_data* input = new dfs_data;
            //Copy data into new struct to pass to new invocation of threadDFS
            input->motif_index = data->motif_index;
            input->cur_node = g_local_nodes[dest_unique_node];
            //Going to a previously visited node, so motif_node_index and role should match
            //every time
            input->desired_motif_node.motif_node_index = next_edge.second.motif_node_index;
            input->desired_motif_node.role = next_edge.second.role;
            input->first_invocation = false;
            input->motif_edges = data->motif_edges;
            input->visited_nodes = data->visited_nodes;
            input->motif_to_unique = data->motif_to_unique;

            threadDFS((void*)input);
            END_DFS(data)
        }
        //Node is in another rank
        else{

        }
    }

    cout << "going to new node" << endl;
    //Know that edge goes to a new node, so iterate over all neighbors
    for(set<int>::iterator neigh = data->cur_node.neighbors.begin();
            neigh != data->cur_node.neighbors.end(); neigh++){
        //Don't need to try any visited nodes
        if(data->visited_nodes.find(*neigh) != data->visited_nodes.end()){
            continue;
        }

        //Determine if the neighbor is in this rank or not
        cout << "before search" << endl;
        rank_with_node = getRankForNode(*neigh);
        cout << "after search" << endl;
        //Node is local, so simply run from current thread
        if(rank_with_node == mpi_myrank){
            cout << "new node" << endl;
            struct dfs_data* input = new dfs_data;
            input->motif_index = data->motif_index;
            input->cur_node = g_local_nodes[*neigh];
            input->desired_motif_node.motif_node_index = next_edge.second.motif_node_index;
            input->desired_motif_node.role = next_edge.second.role;
            input->first_invocation = false;
            input->motif_edges = data->motif_edges;
            input->visited_nodes = data->visited_nodes;
            input->motif_to_unique = data->motif_to_unique;

            threadDFS((void*)input);
        }
        //Node is in another rank
        else{

        }
    }

    cout << "got to end" << endl;
    END_DFS(data)
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
        data->desired_motif_node = data->motif_edges.front().first;

        thpool_add_work(g_local_threads, threadDFS, (void*)data);
        //Test code
        //sleep((mpi_myrank+1) * 1);
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
    //Single rank test
    struct graph_node v;
    v.role = 0;
    v.unique_id = 0;
    g_local_nodes[0] = v;
    v.unique_id = 1;
    g_local_nodes[1] = v;
    v.unique_id = 2;
    v.neighbors.insert(0);
    v.neighbors.insert(1);
    g_local_nodes[2] = v;
    v.neighbors.clear();
    v.unique_id = 3;
    v.neighbors.insert(2);
    v.neighbors.insert(4);
    v.neighbors.insert(7);
    g_local_nodes[3] = v;
    v.neighbors.clear();
    v.unique_id = 4;
    v.neighbors.insert(5);
    v.neighbors.insert(6);
    g_local_nodes[4] = v;
    v.neighbors.clear();
    v.unique_id = 5;
    g_local_nodes[5] = v;
    v.unique_id = 6;
    g_local_nodes[6] = v;
    v.unique_id = 7;
    v.neighbors.insert(8);
    v.neighbors.insert(9);
    g_local_nodes[7] = v;
    v.neighbors.clear();
    v.unique_id = 8;
    v.neighbors.insert(9);
    g_local_nodes[8] = v;
    v.neighbors.clear();
    v.unique_id = 9;
    v.neighbors.insert(10);
    v.neighbors.insert(11);
    g_local_nodes[9] = v;
    v.neighbors.clear();
    v.unique_id = 10;
    v.neighbors.insert(11);
    g_local_nodes[10] = v;
    v.neighbors.clear();
    v.unique_id = 11;
    v.neighbors.insert(3);
    v.neighbors.insert(12);
    g_local_nodes[11] = v;
    v.neighbors.clear();
    v.unique_id = 12;
    v.neighbors.insert(13);
    v.neighbors.insert(15);
    g_local_nodes[12] = v;
    v.neighbors.clear();
    v.unique_id = 13;
    v.neighbors.insert(14);
    g_local_nodes[13] = v;
    v.neighbors.clear();
    v.unique_id = 14;
    g_local_nodes[14] = v;
    v.unique_id = 15;
    v.neighbors.insert(14);
    g_local_nodes[15] = v;

    g_vtxdist.push_back(0);
    g_vtxdist.push_back(16);
    //End test code

    //Test code: create dummy motifs
    //Single rank test
    list<pair<struct motif_node, struct motif_node> > m;
    struct motif_node a;
    a.role = 0;
    struct motif_node b;
    b.role = 0;
    //Motif 0
    a.motif_node_index = 0;
    b.motif_node_index = 1;
    m.push_back(make_pair(a, b));
    b.motif_node_index = 2;
    m.push_back(make_pair(a, b));
    a.motif_node_index = 1;
    m.push_back(make_pair(a, b));
    g_motifs.push_back(m);
    m.clear();
    //Motif 1
    for(int i = 0; i < 4; i++){
        a.motif_node_index = i;
        b.motif_node_index = (i + 1) % 4;
        m.push_back(make_pair(a, b));
    }
    g_motifs.push_back(m);
    m.clear();
    //Motif 2
    for(int i = 0; i < 5; i++){
        a.motif_node_index = i;
        b.motif_node_index = (i + 1) % 5;
        m.push_back(make_pair(a, b));
    }
    g_motifs.push_back(m);
    m.clear();
    //Motif 3
    a.motif_node_index = 0;
    b.motif_node_index = 1;
    m.push_back(make_pair(a, b));
    a.motif_node_index = 1;
    b.motif_node_index = 2;
    m.push_back(make_pair(a, b));
    b.motif_node_index = 3;
    m.push_back(make_pair(a, b));
    g_motifs.push_back(m);
    m.clear();
    //Motif 4
    a.motif_node_index = 0;
    b.motif_node_index = 1;
    m.push_back(make_pair(a, b));
    a.motif_node_index = 1;
    b.motif_node_index = 2;
    m.push_back(make_pair(a, b));
    b.motif_node_index = 3;
    m.push_back(make_pair(a, b));
    a.motif_node_index = 2;
    b.motif_node_index = 4;
    m.push_back(make_pair(a, b));
    a.motif_node_index = 3;
    m.push_back(make_pair(a, b));
    g_motifs.push_back(m);
    m.clear();
    //Motif 5
    a.motif_node_index = 0;
    b.motif_node_index = 1;
    m.push_back(make_pair(a, b));
    m.push_back(make_pair(b, a));
    g_motifs.push_back(m);

    for(unsigned int i = 0; i < g_motifs.size(); i++){
        g_motif_counts.push_back(0);
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
    for(unsigned int i = 0; i < g_motif_counts.size(); i++){
        cout << "Count for motif " << i << ": " << g_motif_counts[i] << endl;
    }
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