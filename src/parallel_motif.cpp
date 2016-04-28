extern "C" {
    #include "C-Thread-Pool/thpool.h"
}
#include <iostream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

void* test_thread(void* id){
    sleep(*(int*)id);
    cout << "Hello from thread " << *(int*)id << " with actual id " << pthread_self() << endl;
    delete (int*)id;
    return NULL;
}

int main(int argc, char* argv[]){
    (void)argc; (void)argv;
    cout << "Test" << endl;
    threadpool tp = thpool_init(4);
    for(int i = 0; i < 12; i++){
        int* id = new int;
        *id = i;
        thpool_add_work(tp, test_thread, (void*)id);
        
    }
    thpool_wait(tp);
    thpool_destroy(tp);
    return 0;
}