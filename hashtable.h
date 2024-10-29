/** \CSE412 Assignment_1
 *  \author Jaeyoung Yun
 *  \sinsunby@unist.ac.kr
 *  \modified Kyu Yeun Kim
 *  \kyuyeunk@unist.ac.kr
 */

// Simply means H_table_t is pointer to 'struct hashtable'
typedef struct hashtable *H_table_t; 

int Hash_Init(H_table_t, int numbuckets);

void Hash_Destroy(H_table_t, int state);

int Hash_Insert(H_table_t, int element);

int Hash_Lookup(H_table_t, int element);

int Hash_Delete(H_table_t, int element);
