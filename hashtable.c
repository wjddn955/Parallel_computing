#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "hashtable.h"

enum operation{ //enumeration of hashtable operation types
    INSERT = 1,	LOOKUP = 2, DELETE = 3
}; 

struct element{ // hashtable element
    struct element *next;
    int key;
};

struct bucket_t{ //hashtable bucket
    struct element *chain;
    pthread_mutex_t mutex; 
};

struct hashtable{ //hashtable
    int num_bucket;
    struct bucket_t *bucket;
};

struct arg_struct{ //passing argument struct
    H_table_t hash_table;
    int thr_id;
};

struct data_t{ //specified data type
    int op; // hash table operation type
    int key; // key value
};

// a pointer variable pointing to the sequence of
// hash table operations performed by multiple threads
struct data_t *data;

int num_op; // total number of operations to be performed
int nthr; // number of threads

int hashtable_init(H_table_t hash_table, int numbuckets){

    //TODO: Initialize the hash table with specified number of buckets.
    hash_table->num_bucket = numbuckets; // Set the num of buckets
    hash_table->bucket = (struct bucket_t *) malloc(numbuckets * sizeof(struct bucket_t)); // Set the bucket to point the memory block(allocated by malloc)

    if (hash_table->bucket == NULL) return -1; // If allocating bucket failed, return -1

    for (int i = 0; i < numbuckets; i++) { // Traverse through buckets and initiate their chain, next, key, initiate mutex
        hash_table->bucket[i].chain = (struct element *) malloc(sizeof (struct element));
        hash_table->bucket[i].chain->next = NULL;
        hash_table->bucket[i].chain->key = -1;
        pthread_mutex_init(&(hash_table->bucket[i].mutex), NULL);
    }

    return 1;
}

int hashtable_insert(H_table_t hash_table, int key){

    //TODO: Insert a key value to the hash table. 
    //If the key value has been successfully added,
    //return 0; otherwise, return -1
    //duplicates are allowed.

    int bucket_idx = key % hash_table->num_bucket; // Calculate bucket index
    struct element *new_element = (struct element *) malloc(sizeof(struct element)); // allocate memory for new_element
    new_element->key = key;

    pthread_mutex_lock(&(hash_table->bucket[bucket_idx].mutex)); // Acquire lock to modify hashtable
    struct element *current = hash_table->bucket[bucket_idx].chain;
    if (current->key == -1) { // Check if chain is empty
        current->key = key;
        free(new_element);
        pthread_mutex_unlock(&(hash_table->bucket[bucket_idx].mutex));
        
        return 0;
    }

    if (new_element == NULL) { // Allocating memory for new_element failed
        free(new_element);
        pthread_mutex_unlock(&(hash_table->bucket[bucket_idx].mutex));
        return -1;
    }
    
    new_element->next = current; // Position new_element before current (Put new element at chain position)
    hash_table->bucket[bucket_idx].chain = new_element;
    pthread_mutex_unlock(&(hash_table->bucket[bucket_idx].mutex));

    return 0;
}

int hashtable_lookup(H_table_t hash_table, int key){

    //TODO: Lookup the key. If it is in the hash table, 
    //return 0; otherwise, return -1

    int ret = -1;
    int bucket_idx = key % hash_table->num_bucket;

    pthread_mutex_lock(&(hash_table->bucket[bucket_idx].mutex)); // Although lookup does not modify hashtable, it loads it so, acquire lock
    struct element *current = hash_table->bucket[bucket_idx].chain;
    while (current != NULL) { // Traverse inside bucket, keep it until the end because duplicate keys are allowed
        if (current->key == key) {
            ret = 0;
        }
        current = current->next;
    } 
    
    pthread_mutex_unlock(&(hash_table->bucket[bucket_idx].mutex));
    return ret;

}

int hashtable_delete(H_table_t hash_table, int key){

    //TODO: Delete all elements with the key.
    //If it is in the hash table,
    //return 0; otherwise, return -1
    int bucket_idx = key % hash_table->num_bucket;
    int ret = -1;
    struct element *prev = NULL;
    
    pthread_mutex_lock(&(hash_table->bucket[bucket_idx].mutex));
    struct element *current = hash_table->bucket[bucket_idx].chain;
    // Traverse until the end, because duplicate keys are allowed
    while (current != NULL) { // Traverse inside bucket and modify the link between elements after delete target element
        if (current->key == key) {
            if (prev == NULL) {
                hash_table->bucket[bucket_idx].chain = current->next;
            }
            else {
                prev->next = current->next;
            }
            struct element *to_delete = current;
            current = current->next;
            free(to_delete);
            ret = 0;
        }
        else {
            prev = current;
            current = current->next;
        }
    }
    pthread_mutex_unlock(&(hash_table->bucket[bucket_idx].mutex));
    
    return ret;
}

void *test(void *arguments){ // pthread Function
    int i, rc, thr_id;
    int start, end;
    H_table_t hash_table;
    //TODO: Assign each variable with a correct value
    //use arguments to get hash table structure and thr_id

    struct arg_struct *args = (struct args_struct *)arguments;
    hash_table = args->hash_table;
    thr_id = args->thr_id;
    
    // Evenly divide workload between threads
    start = (num_op/nthr) * thr_id;
    if (thr_id == nthr -1) {
        end = num_op;
    } else {
        end = start + (num_op / nthr);
    }

    for(i= start; i< end; i++){
        if(data[i].op == INSERT)
            rc=hashtable_insert(hash_table, data[i].key);
        else if(data[i].op == LOOKUP)
            rc=hashtable_lookup(hash_table, data[i].key);
        else if(data[i].op == DELETE)
            rc=hashtable_delete(hash_table, data[i].key);
    }

    free(args);
    return 0;
}


int main(int argc, char** argv){
    struct timespec time_info;
    const int hash_num_bucket = atoi(argv[1]);
    num_op = atoi(argv[2]);
    nthr = atoi(argv[3]);   
    const int i_ratio = atoi(argv[4]);
    const int d_ratio = atoi(argv[5]);
    int i, thr_id, status, result=0;
    int64_t start_time, end_time;
    H_table_t hash_table;

    pthread_t *p_thread = malloc(nthr * sizeof(pthread_t));
    if (p_thread == NULL) {
        perror("malloc p_thread");
        exit(EXIT_FAILURE);
    }
    
    struct arg_struct *args = malloc(nthr * sizeof(struct arg_struct));
    if (args == NULL) {
        perror("malloc args");
        exit(EXIT_FAILURE);
    } 

    // 1. Data generation phase
    data = malloc(num_op * sizeof(struct data_t));	
    if (data == NULL) {
        perror("malloc data");
        exit(EXIT_FAILURE);
    } 
    
    srand(22);
    for(i=0; i<num_op; i++){ //make a data set
        int r = rand()%65536;

        data[i].key = r;
        r = rand()%100;
        if(r < i_ratio)
            data[i].op = INSERT;
        else if(r < i_ratio + d_ratio)
            data[i].op = DELETE;
        else
            data[i].op = LOOKUP;
    }
    printf("data success\n");

    hash_table = (H_table_t)malloc(sizeof(struct hashtable));
    if (hash_table == NULL) {
        perror("malloc hash_table");
        exit(EXIT_FAILURE);
    } 
    
    hashtable_init(hash_table, hash_num_bucket); //Initialize the hash table
    printf("init success\n");
    
    // 2. Performance test
    if(	clock_gettime( CLOCK_REALTIME, &time_info ) == -1 ){ //record start point
        perror( "clock gettime" );
        exit( EXIT_FAILURE );
    }
    
    start_time = ( (int64_t) time_info.tv_sec * 1000000000 + (int64_t) time_info.tv_nsec );
    
    for(i=0; i<nthr; i++){ //make threads and pass the arguments
        //TODO: create thread by pthread_create()
        // use arg_struct to pass the multiple arguments
        // (thread_id, hash_table)

        // Make arg to pass to threads
        struct arg_struct *arg = malloc(sizeof(struct arg_struct));
        if (arg == NULL) {
            perror("malloc arg");
            exit(EXIT_FAILURE);
        }
        arg->hash_table = hash_table;
        arg->thr_id = i;
        status = pthread_create(&p_thread[i], NULL, test, (void *)arg);
        if (status != 0) {
            exit(EXIT_FAILURE);
        }
    }
    for(i=0; i<nthr; i++) {
        //TODO: use pthread_join to wait till the thread terminate

        // Use pthread_join
        status = pthread_join(p_thread[i], NULL);
        if (status != 0) exit(EXIT_FAILURE);
    }
    if(	clock_gettime( CLOCK_REALTIME, &time_info ) == -1 ){ //record end point
            perror( "clock gettime" );
            exit( EXIT_FAILURE );
    }
    end_time = ( (int64_t) time_info.tv_sec * 1000000000 + (int64_t) time_info.tv_nsec );

    printf("Running Time: %lf s\n", (double)((long)end_time-(long)start_time)/1000000000.0);

    // 3. Correctness test
    // DO NOT MODIFY THE CODE BELOW

    for(i=0; i<num_op; i++){ //simple test for the Hash table
        bool inserted=false;
        bool deleted=false;
        int j;
        for(j=0; j<num_op; j++){
            if(data[j].key==data[i].key){
                if(data[j].op==INSERT){
                    inserted=true;
                } 
                else if(data[j].op==DELETE){
                    deleted=true;
                }    
            }    
        }
        
        int lookup_result=hashtable_lookup(hash_table, data[i].key);
        if((inserted==true && deleted==false && lookup_result==-1) || 
           (inserted==false && lookup_result==0) || 
           lookup_result==1){
            result=-1;
            printf("error 1:%d\n",i);    
        }
    }   

    //result of test
    if(result==-1)
        printf("Correctness: Incorrect\n");
    else
        printf("Correctness: Correct\n");
    return 0;
}
