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
    #include <semaphore.h>
}

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <list>
#include <string>
#include <mpi.h> //MPI apparently does the extern stuff in its header file

//MPI tags
#define MPIT_RANK_DONE 0
#define MPIT_REQUEST_NODE 1
#define MPIT_RESPONSE_NODE 2
#define MPIT_EDGE_START 3
#define MPIT_EDGE_END 4

#define CLOCKRATE 1600000000

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
    map<int, struct graph_node> remote_cache;
};

struct node_request { //Struct sent over MPI when one rank requests node data from another
    pthread_t thread_id;
    int node;
};

struct node_response { //struct sent over MPI responding to a request for node data
    pthread_t thread_id;
    int unique_id;
    int role;
    //neighbors data is sent in the space in the buffer after this struct
};

/***************************************************************************/
/* Global Variables ********************************************************/
/***************************************************************************/
string dataFile;
string motifFile;
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
pthread_mutex_t* m_shared_buffers;
pthread_mutex_t* m_locked_threads;
map<pthread_t, sem_t*> g_locked_threads;
map<pthread_t, struct graph_node*> g_shared_buffers; //Used to pass received nodes to threads
                                                     //Thread is responsible for freeing memory

//Timing
unsigned long long g_gen_time_start = 0;
unsigned long long g_gen_time_end = 0;
unsigned long long g_comp_time_start = 0;
unsigned long long g_comp_time_end = 0;

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
    
        cerr << "Requested to find a node not in the graph " << id<<endl;
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
            //cout << "Found rank " << m << " for node " << id << endl;
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
/* getNodeFromRank *********************************************************/
/***************************************************************************/
//Sends a request for the given node located in the given rank and stores it in
//the given cache
void getNodeFromRank(int node_id, int rank, map<int, struct graph_node> &cache){
    pthread_t thread_id = pthread_self();

    //Create a semaphore for blocking until a response is received
    sem_t* sem = new sem_t;
    int rc = sem_init(sem, 0, 0); //Initialize a shared semaphore with initial value 0
    if(rc != 0){
        perror("Semaphore initialization failed");
        CLEAN_EXIT
    }
    pthread_mutex_lock(m_locked_threads);
    g_locked_threads[thread_id] = sem;
    pthread_mutex_unlock(m_locked_threads);

    //Send the request
    MPI_Request req = MPI_REQUEST_NULL;
    struct node_request node_req;
    node_req.thread_id = thread_id;
    node_req.node = node_id;
    pthread_mutex_lock(m_mpi_lock);
    MPI_Isend(&node_req, sizeof(node_req), MPI_BYTE, rank, MPIT_REQUEST_NODE, MPI_COMM_WORLD, &req);
    pthread_mutex_unlock(m_mpi_lock);

    //Block until main thread unlocks this thread
    sem_wait(sem);

    //Copy node into cache
    pthread_mutex_lock(m_shared_buffers);
    cache[node_id] = *g_shared_buffers[thread_id];
    delete g_shared_buffers[thread_id];
    pthread_mutex_unlock(m_shared_buffers);

    rc = sem_destroy(sem);
    if(rc != 0){
        perror("Semaphore destruction failed");
        CLEAN_EXIT
    }
    delete sem;
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
    //Need a fresh copy of motif edges
    data->motif_edges = g_motifs[data->motif_index];

    //Go through each of the visited nodes
    for(map<int, int>::iterator m_to_u_iter = data->motif_to_unique.begin();
            m_to_u_iter != data->motif_to_unique.end(); m_to_u_iter++){
        int cur_node_m_id = m_to_u_iter->first;
        struct graph_node cur_node;
        map<int, struct graph_node>::iterator node_iter = g_local_nodes.find(m_to_u_iter->second);
        //Remote node
        if(node_iter == g_local_nodes.end()){
            cur_node = data->remote_cache[m_to_u_iter->second];
        }
        //Local node
        else{
            cur_node = node_iter->second;
        }
        set<int> neighbors_in_motif;

        //Get all the nodes in the motif that this node has a motif edge to
        for(list<pair<struct motif_node, struct motif_node> >::iterator edge_iter = data->motif_edges.begin();
                edge_iter != data->motif_edges.end(); edge_iter++){
            if(edge_iter->first.motif_node_index != cur_node_m_id){
                continue;
            }
            neighbors_in_motif.insert(data->motif_to_unique[edge_iter->second.motif_node_index]);
        }

        //See if there are any non-motif edges between this node and motif neighbors
        set<int>::iterator neighbor_iter;
        for(set<int>::iterator visited_iter = data->visited_nodes.begin();
                visited_iter != data->visited_nodes.end(); visited_iter++){

            neighbor_iter = cur_node.neighbors.find(*visited_iter);
            //The visited node is not a neighbor of the node in question
            if(neighbor_iter == cur_node.neighbors.end()){
                continue;
            }

            //We know that there exists an edge between the node in question and the visited node
            //If there isn't an edge between the two in the motif, this is not a valid motif
            neighbor_iter = neighbors_in_motif.find(*visited_iter);
            if(neighbor_iter == neighbors_in_motif.end()){
                return false;
            }
        }
    }
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
            get remote node data and call threadDFS
            return
    (next edge starts at the current node)
    Pop the next edge and store it
    if edge goes to a previously visited node: (Optimization for closing a cycle, as we already know which node we need to go to)
        if that node is a neighbor:
            if neighbor is local:
                call theadDFS on that node and wait
                return
            else:
                get remote node data and call threadDFS
                return
    (Know that the node we're looking for has not yet been visited)
    for each neighbor:
        if neighbor has been visited:
            continue
        if neighbor is local:
            call threadDFS on that node and wait
            
        else:
            get remote node data and call threadDFS
            remove remote node from cache
            */
void* threadDFS(void* d){
    struct dfs_data* data = (struct dfs_data*)d;

    //Test code
    /*cout << "Rank " << mpi_myrank << " DFS run on motif " << data->motif_index << " graph node " << data->cur_node.unique_id;
    cout << " desired motif node role " << data->desired_motif_node.motif_node_index << " motif edges ";
    for(list<pair<struct motif_node, struct motif_node> >::iterator i = data->motif_edges.begin();
            i != data->motif_edges.end(); i++){
        cout << "([" << i->first.motif_node_index << ", ";
        cout << i->first.role << "], [" << i->second.motif_node_index;
        cout << ", " << i->second.role << "]) ";
    }
    cout << endl;*/
    //End test code

    //Ensure we're at a valid node for the motif
    //Current node and desired motif node role mismatch
    if(data->cur_node.role != data->desired_motif_node.role){
        //cout << "mismatched role" << endl;
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

    //cout << "valid node" << endl;
    //Now know that the node we're at is valid
    //Check if this was the last node needed to complete motif
    if(data->motif_edges.empty()){
        if(isValidMotif(data)){
            pthread_mutex_lock(m_counter_lock);
            g_motif_counts[data->motif_index] += 1;
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

    //cout << "valid motif representation" << endl;
    //next_node is a unique node identifier
    int source_unique_node = it->second;
    int rank_with_node;
    //Need to hop to previously visited node before continuing
    if(source_unique_node != data->cur_node.unique_id){
        //Determine if the visited node is in this rank or not
        rank_with_node = getRankForNode(source_unique_node);
        //Node is local, so simply jump
        if(rank_with_node == mpi_myrank){
            //cout << "jumping" << endl;
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
            input->remote_cache = data->remote_cache;

            threadDFS((void*)input);
            END_DFS(data)
        }
        //Node is in another rank
        else{
            //Should be in the cache but make sure
            map<int, struct graph_node>::iterator cache_it = data->remote_cache.find(source_unique_node);
            if(cache_it == data->remote_cache.end()){
                cout << "Cache miss that shouldn't happen 1" << endl;
                getNodeFromRank(source_unique_node, rank_with_node, data->remote_cache);
            }

            struct dfs_data* input = new dfs_data;
            input->motif_index = data->motif_index;
            input->cur_node = data->remote_cache[source_unique_node];
            //Jumping, so motif_node_index should be the same as the previously visited node
            input->desired_motif_node.motif_node_index = it->first;
            //Role should be the same as the previously visited node
            input->desired_motif_node.role = input->cur_node.role;
            input->first_invocation = false;
            input->motif_edges = data->motif_edges;
            input->visited_nodes = data->visited_nodes;
            input->motif_to_unique = data->motif_to_unique;
            input->remote_cache = data->remote_cache;
            threadDFS((void*)input);
            END_DFS(data);
        }
    }

    //cout << "continues from current node" << endl;
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
            //cout << "previously visited node not neighbor" << endl;
            END_DFS(data)
        }
        //Determine if the visited node is in this rank or not
        rank_with_node = getRankForNode(dest_unique_node);
        //cout << "Was given rank " << rank_with_node << " for node " << dest_unique_node << endl;
        //Node is local, so simply call from same thread
        if(rank_with_node == mpi_myrank){
            //cout << "previously visited" << endl;
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
            input->remote_cache = data->remote_cache;

            threadDFS((void*)input);
            END_DFS(data)
        }
        //Node is in another rank
        else{
            //Should be in cache, but make sure
            map<int, struct graph_node>::iterator cache_it = data->remote_cache.find(dest_unique_node);
            if(cache_it == data->remote_cache.end()){
                cout << "Cache miss that shouldn't happen 2 for node " << dest_unique_node << " in rank " << mpi_myrank << endl;
                getNodeFromRank(dest_unique_node, rank_with_node, data->remote_cache);
            }

            struct dfs_data* input = new dfs_data;
            input->motif_index = data->motif_index;
            input->cur_node = data->remote_cache[dest_unique_node];
            //Going to previously visited node, so index and role should match every time
            input->desired_motif_node.motif_node_index = next_edge.second.motif_node_index;
            input->desired_motif_node.role = next_edge.second.role;
            input->first_invocation = false;
            input->motif_edges = data->motif_edges;
            input->visited_nodes = data->visited_nodes;
            input->motif_to_unique = data->motif_to_unique;
            input->remote_cache = data->remote_cache;

            threadDFS((void*)input);
            END_DFS(data);
        }
    }

    //cout << "going to new node" << endl;
    //Know that edge goes to a new node, so iterate over all neighbors
    for(set<int>::iterator neigh = data->cur_node.neighbors.begin();
            neigh != data->cur_node.neighbors.end(); neigh++){
        //Don't need to try any visited nodes
        if(data->visited_nodes.find(*neigh) != data->visited_nodes.end()){
            continue;
        }

        //Determine if the neighbor is in this rank or not
        //cout << "before search" << endl;
        rank_with_node = getRankForNode(*neigh);
        //cout << "after search" << endl;
        //Node is local, so simply run from current thread
        if(rank_with_node == mpi_myrank){
            //cout << "new node" << endl;
            struct dfs_data* input = new dfs_data;
            input->motif_index = data->motif_index;
            input->cur_node = g_local_nodes[*neigh];
            input->desired_motif_node.motif_node_index = next_edge.second.motif_node_index;
            input->desired_motif_node.role = next_edge.second.role;
            input->first_invocation = false;
            input->motif_edges = data->motif_edges;
            input->visited_nodes = data->visited_nodes;
            input->motif_to_unique = data->motif_to_unique;
            input->remote_cache = data->remote_cache;

            threadDFS((void*)input);
        }
        //Node is in another rank
        else{
            //Going to a new node, so should not be in cache
            map<int, struct graph_node>::iterator cache_it = data->remote_cache.find(*neigh);
            if(cache_it != data->remote_cache.end()){
                cout << "Cache hit that shouldn't happen" << endl;
            }

            getNodeFromRank(*neigh, rank_with_node, data->remote_cache);

            struct dfs_data* input = new dfs_data;
            input->motif_index = data->motif_index;
            input->cur_node = data->remote_cache[*neigh];
            input->desired_motif_node.motif_node_index = next_edge.second.motif_node_index;
            input->desired_motif_node.role = next_edge.second.role;
            input->first_invocation = false;
            input->motif_edges = data->motif_edges;
            input->visited_nodes = data->visited_nodes;
            input->motif_to_unique = data->motif_to_unique;
            input->remote_cache = data->remote_cache;

            threadDFS((void*)input);

            //Remove node from cache so cache doesn't balloon during loop
            data->remote_cache.erase(*neigh);
        }
    }

    //cout << "got to end" << endl;
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

void genTestData();
void printStartInfo();

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
    dataFile=string(argv[1]);
    motifFile=string(argv[2]);
    //Initialize thread pools
    g_local_threads = thpool_init(num_local_threads);

    //Initialize synchronization objects
    m_counter_lock = new pthread_mutex_t;
    rc = pthread_mutex_init(m_counter_lock, NULL);
    CHECK_MUTEX_INIT(rc)
    m_mpi_lock = new pthread_mutex_t;
    rc = pthread_mutex_init(m_mpi_lock, NULL);
    CHECK_MUTEX_INIT(rc)
    m_shared_buffers = new pthread_mutex_t;
    rc = pthread_mutex_init(m_shared_buffers, NULL);
    CHECK_MUTEX_INIT(rc)
    m_locked_threads = new pthread_mutex_t;
    rc = pthread_mutex_init(m_locked_threads, NULL);
    CHECK_MUTEX_INIT(rc)

    if(mpi_myrank == 0){
        g_gen_time_start = GetTimeBase();
    }
    genTestData();
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_myrank == 0){
        g_gen_time_end = GetTimeBase();
    }
    sleep(mpi_myrank * 2);
    printStartInfo();
    MPI_Barrier(MPI_COMM_WORLD);

    if(mpi_myrank == 0){
        g_comp_time_start = GetTimeBase();
    }

    //Search for all motifs
    MPI_Status status;
    int flag;
    int* motif_index = new int;
    for(unsigned int i = 0; i < g_motifs.size(); i++){
    //cout<<"at motif 1"<<endl;
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
                rc = MPI_Recv(NULL, 0, MPI_INT, status.MPI_SOURCE, MPIT_RANK_DONE, 
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                pthread_mutex_unlock(m_mpi_lock);
                if(rc != MPI_SUCCESS){
                    cerr << "Failed to clear rank done message with error" << rc << endl;
                    CLEAN_EXIT
                }
                g_ranks_done++;
                if(g_ranks_done == mpi_commsize){
                    break; //All ranks are done, so break out of infinite loop
                }
            }
            //Received a request for node data
            else if(status.MPI_TAG == MPIT_REQUEST_NODE){
                //Respond to request
                //Grab the incoming message
                struct node_request node_req;
                pthread_mutex_lock(m_mpi_lock);
                rc = MPI_Recv(&node_req, sizeof(node_req), MPI_BYTE, status.MPI_SOURCE, MPIT_REQUEST_NODE,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                pthread_mutex_unlock(m_mpi_lock);

                map<int, struct graph_node>::iterator iter = g_local_nodes.find(node_req.node);
                if(iter == g_local_nodes.end()){
                    cerr << "Rank " << status.MPI_SOURCE << " requested non-existent node ";
                    cerr << node_req.node << " from rank " << mpi_myrank << endl;
                    CLEAN_EXIT
                }

                //Prepare struct for sending
                struct node_response node_resp;
                node_resp.thread_id = node_req.thread_id;
                node_resp.unique_id = iter->second.unique_id;
                node_resp.role = iter->second.role;

                //Create an exact size buffer
                int buf_size = sizeof(node_resp) + sizeof(int) * iter->second.neighbors.size();
                char send_buffer[buf_size];

                //Copy data into buffer
                memcpy(send_buffer, &node_resp, sizeof(node_resp));
                set<int>::iterator neigh_iter = iter->second.neighbors.begin();
                int elem;
                for(int i = 0; neigh_iter != iter->second.neighbors.end(); i++, neigh_iter++){
                    elem = *neigh_iter;
                    memcpy(&send_buffer[sizeof(node_resp) + i * sizeof(elem)], &elem, sizeof(elem));
                }

                //Send buffer
                pthread_mutex_lock(m_mpi_lock);
                rc = MPI_Send(send_buffer, buf_size, MPI_BYTE, status.MPI_SOURCE, 
                        MPIT_RESPONSE_NODE, MPI_COMM_WORLD);
                pthread_mutex_unlock(m_mpi_lock);
                if(rc != MPI_SUCCESS){
                    cerr << "Failed to send node response with error code " << rc << endl;
                    CLEAN_EXIT
                }
                /*cout << "Sent node " << node_resp.unique_id << " with neighbors ";
                for(neigh_iter = iter->second.neighbors.begin(); neigh_iter != iter->second.neighbors.end(); neigh_iter++){
                    cout << *neigh_iter << " ";
                }
                cout << endl;*/
            }
            //Received a response to a request for node data
            else if(status.MPI_TAG == MPIT_RESPONSE_NODE){
                //Copy into struct
                //Get the size of the incoming message
                int num_bytes;
                pthread_mutex_lock(m_mpi_lock);
                rc = MPI_Get_count(&status, MPI_BYTE, &num_bytes);
                pthread_mutex_unlock(m_mpi_lock);
                if(rc != MPI_SUCCESS){
                    cerr << "Failed to get number of bytes with error code " << rc << endl;
                    CLEAN_EXIT
                }

                //cout << "Received num_bytes " << num_bytes << " which should have " << (num_bytes - sizeof(struct node_response)) / sizeof(int) << " neighbors" << endl;

                //Receive data
                char recv_buffer[num_bytes];
                pthread_mutex_lock(m_mpi_lock);
                rc = MPI_Recv(recv_buffer, num_bytes, MPI_BYTE, status.MPI_SOURCE,
                        MPIT_RESPONSE_NODE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                pthread_mutex_unlock(m_mpi_lock);
                if(rc != MPI_SUCCESS){
                    cerr << "Failed to receive node response with error code " << rc << endl;
                    CLEAN_EXIT
                }

                //Interpret data
                struct node_response node_resp;
                struct graph_node* g = new struct graph_node;
                memcpy(&node_resp, recv_buffer, sizeof(node_resp));
                g->unique_id = node_resp.unique_id;
                g->role = node_resp.role;
                int temp;
                for(int i = sizeof(node_resp); i < num_bytes; i += sizeof(int)){
                    memcpy(&temp, &recv_buffer[i], sizeof(int));
                    //cout << "Copied element " << temp << endl;
                    g->neighbors.insert(temp);
                }

                /*cout << "Received node " << node_resp.unique_id << " with neighbors " ;
                for(set<int>::iterator iter = g->neighbors.begin(); iter != g->neighbors.end(); iter++){
                    cout << *iter << " ";
                }
                cout << endl;*/

                //Store pointer for thread and unlock thread
                pthread_mutex_lock(m_shared_buffers);
                g_shared_buffers[node_resp.thread_id] = g;
                pthread_mutex_unlock(m_shared_buffers);

                pthread_mutex_lock(m_locked_threads);
                rc = sem_post(g_locked_threads[node_resp.thread_id]);
                pthread_mutex_unlock(m_locked_threads);
                if(rc != 0){
                    perror("Failed to unlock thread");
                    CLEAN_EXIT
                }
            }
            else{
                cout << "Received unknown tag " << status.MPI_TAG << endl;
                CLEAN_EXIT
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        //Reduce all counts to rank 0
        int global_sum;
        //cout << "Rank " << mpi_myrank << " found " << g_motif_counts[i] << " for motif " << i << endl;
        rc = MPI_Reduce(&g_motif_counts[i], &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(rc != MPI_SUCCESS){
            cerr << "Failed to reduce with error code " << rc << endl;
            CLEAN_EXIT
        }
        if(mpi_myrank == 0){
            g_motif_counts[i] = global_sum;
        }
    }
    delete motif_index;
    CLEANUP_TP

    if(mpi_myrank == 0){
        g_comp_time_end = GetTimeBase();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    rc = MPI_Finalize();
    if(mpi_myrank == 0){
        for(unsigned int i = 0; i < g_motif_counts.size(); i++){
            cout << "Count for motif " << i << ": " << g_motif_counts[i] << endl;
        }
        cout << "Graph distribution time took " << (g_gen_time_end - g_gen_time_start) / ((double)CLOCKRATE) << " seconds" << endl;
        cout << "Computation time took " << (g_comp_time_end - g_comp_time_start) / ((double)CLOCKRATE) << " seconds" << endl;
    }    
    return 0;
}

/***************************************************************************/
/* genTestData *************************************************************/
/***************************************************************************/

void genTestData(){
    int totalNodes, nodesPerRank;
    if(mpi_myrank==0){
       
        std::fstream myfile(dataFile.c_str(), std::ios_base::in);
        //cout<<"file loc"<<dataFile<<endl;
        
        //read in totalnumber of nodes
        myfile>>totalNodes;
        
        //Let every rank know total number of nodes
        MPI_Bcast(&totalNodes,1,MPI_INT,0, MPI_COMM_WORLD);
        nodesPerRank=totalNodes/mpi_commsize;

        //Generate vtxdist
        for(int i = 0; i <= totalNodes; i += nodesPerRank){
            if(i+nodesPerRank>totalNodes){
                g_vtxdist.push_back(totalNodes);
            }
            else{
                g_vtxdist.push_back(i);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        int desRank;
        int edgeVal[2]; 
        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Status status;
        std::map<int,struct graph_node>::iterator it;
        struct graph_node v;
        v.role = 0;
    
        while (myfile >> edgeVal[0] >> edgeVal[1]){
        
            //Add the source node
            desRank = getRankForNode(edgeVal[0]);
            if(desRank != 0){
                MPI_Isend(&edgeVal, 2, MPI_INT, desRank, MPIT_EDGE_START, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
            }
            else{
                it = g_local_nodes.find(edgeVal[0]);
                if(it == g_local_nodes.end()){
                    v.neighbors.clear();
                    v.unique_id = edgeVal[0];
                    v.neighbors.insert(edgeVal[1]);
                    g_local_nodes[edgeVal[0]] = v;
                }
                else{
                    it->second.neighbors.insert(edgeVal[1]);
                }
            }

            //Add the destination node. Necessary in case of sink nodes
            desRank = getRankForNode(edgeVal[1]);
            if(desRank != 0){
                MPI_Isend(&edgeVal, 2, MPI_INT, desRank, MPIT_EDGE_END, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
            }
            else{
                it = g_local_nodes.find(edgeVal[1]);
                if(it == g_local_nodes.end()){
                    v.neighbors.clear();
                    v.unique_id = edgeVal[1];
                    g_local_nodes[edgeVal[1]] = v;
                }
            }
        }
    
        myfile.close();
        //indicate that the edges have been read
        edgeVal[0]=-1;
        edgeVal[1]=-1;
        for(int i=1;i<mpi_commsize;i++){
            MPI_Isend(&edgeVal, 2, MPI_INT, i, MPIT_EDGE_START, MPI_COMM_WORLD, &req);
            MPI_Wait(&req, &status);
        }
    }
    else{    
        //recieve total number of nodes and generate vtxdist
        MPI_Bcast(&totalNodes,1,MPI_INT,0, MPI_COMM_WORLD);
    
        nodesPerRank=totalNodes/mpi_commsize;
        
        //Generate vtxdist
        for(int i = 0; i <= totalNodes; i += nodesPerRank){
            if(i+nodesPerRank>totalNodes){
                g_vtxdist.push_back(totalNodes);
            }
            else{
                g_vtxdist.push_back(i);
            }
        }
    
        MPI_Barrier(MPI_COMM_WORLD);

        std::map<int,struct graph_node>::iterator it;
        int edgeVal[2];
        edgeVal[0]=0;
        edgeVal[1]=1; 
        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Status status;
        struct graph_node v;
        v.role = 0;
       
        while(edgeVal[0] != -1){
            MPI_Irecv(&edgeVal, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,&req);
            MPI_Wait(&req, &status);
            if(status.MPI_TAG == MPIT_EDGE_START){
                if(edgeVal[0] != -1){
                    it=g_local_nodes.find(edgeVal[0]);
                    if(it==g_local_nodes.end()){
                        v.neighbors.clear();
                        v.unique_id = edgeVal[0];
                        v.neighbors.insert(edgeVal[1]);
                        g_local_nodes[edgeVal[0]]=v;
                    }
                    else{
                        it->second.neighbors.insert(edgeVal[1]);
                    }
                }
            }
            else{
                it = g_local_nodes.find(edgeVal[1]);
                if(it == g_local_nodes.end()){
                    v.neighbors.clear();
                    v.unique_id = edgeVal[1];
                    g_local_nodes[edgeVal[1]] = v;
                }
            }
        }
    }

    //printf("getmotifs\n");
    list<pair<struct motif_node, struct motif_node> > m;
    struct motif_node a;
    a.role = 0;
    struct motif_node b;
    b.role = 0;
    
    if(mpi_myrank==0){
        std::fstream myfile(motifFile.c_str(), std::ios_base::in);
        //cout<<"file loc "<<motifFile<<endl;
        int edgeVal[2];
            
        while (myfile >> edgeVal[0] >> edgeVal[1]){
            if(edgeVal[0]!=-1){
                //printf("new edge %d %d \n",edgeVal[0],edgeVal[1]);
                a.motif_node_index = edgeVal[0];
                b.motif_node_index = edgeVal[1];
                m.push_back(make_pair(a, b));
                MPI_Bcast(&edgeVal,2,MPI_INT,0, MPI_COMM_WORLD);
            }
            if(edgeVal[0]==-1){
                g_motifs.push_back(m);
                m.clear();
                edgeVal[0]=-1;
                edgeVal[1]=-1;
                //printf("end motif\n");
                MPI_Bcast(&edgeVal,2,MPI_INT,0, MPI_COMM_WORLD);
            }
        }

        g_motifs.push_back(m);
        m.clear();
        //printf("end of file \n");
        edgeVal[0]=-2;
        edgeVal[1]=-2;
        MPI_Bcast(&edgeVal,2,MPI_INT,0, MPI_COMM_WORLD);
        
        myfile.close();
    }
    else{
        int edgeVal[2];
        MPI_Bcast(&edgeVal,2,MPI_INT,0, MPI_COMM_WORLD);

        while(edgeVal[0]!=-2){
            if(edgeVal[0]!=-1){
                //printf("new edge %d %d \n",edgeVal[0],edgeVal[1]);
                a.motif_node_index = edgeVal[0];
                b.motif_node_index = edgeVal[1];
                m.push_back(make_pair(a, b));
            }
            else{
                //printf("new motif\n");
                g_motifs.push_back(m);
                m.clear();
            }
            MPI_Bcast(&edgeVal,2,MPI_INT,0, MPI_COMM_WORLD);
        }

        g_motifs.push_back(m);
    }
    //printf(" rank %d  has %d motifs to llok for \n",mpi_myrank,g_motifs.size());
     for(unsigned int i = 0; i < g_motifs.size(); i++){
        g_motif_counts.push_back(0);
    }
/*
    //Test code: create dummy motifs
    //Single rank test
    
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
    //End test code*/
}

void printStartInfo(){
    printGraph();
    cout << "vtxdist: ";
    for(unsigned int i = 0; i < g_vtxdist.size(); i++){
        cout << g_vtxdist[i] << " ";
    }
    cout << endl;
}
