#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>

int* merge(int *A, int *B, int asize, int bsize);
void mergesort(int *A, int min, int max);

double start_time, stop_time;

int* merge(int *A, int *B, int asize, int bsize) {
	// TODO: fill in the code here to merge the sorted arrays

	// Allocate new array with (asize + bsize) size
	int* new_arr = (int *)malloc((asize + bsize) * sizeof(int));
	int i = 0, j = 0, k = 0;

	// Do sort with assuming A and B are sorted
	while (j < asize && k < bsize) {
        if (A[j] <= B[k]) { new_arr[i++] = A[j++]; } 
		else { new_arr[i++] = B[k++]; }
    }

    // Copy remaining elements from A or B
    while (j < asize) { new_arr[i++] = A[j++]; }
    while (k < bsize) { new_arr[i++] = B[k++]; }

    return new_arr;
}

void mergesort(int *A, int min, int max)
{
	// TODO: fill in the code here to recursive divide the array
	// into two halves, sort and merge them
	
	// Immediately return if mergesort is impossible
	if (min >= max) return;

	// Caldulate mid point (idx)
	int mid = (min + max) / 2;

	// Do mergesort for left half and right half, as a result, the left and right parts got sorted (Not entire)
	mergesort(A, min, mid);
	mergesort(A, mid + 1, max);

	int left_size = mid - min + 1;
    int right_size = max - mid;

    int *left = (int *)malloc(left_size * sizeof(int));
    int *right = (int *)malloc(right_size * sizeof(int));

	for (int i = 0; i < left_size; i++) {
        left[i] = A[min + i];
    }
    for (int i = 0; i < right_size; i++) {
        right[i] = A[mid + 1 + i];
    }

	// Call merge function to merge left and right part
    int *merged = merge(left, right, left_size, right_size);

	// Copy merged data to A
	for (int i = 0; i < left_size + right_size; i++) {
        A[min + i] = merged[i];
    }

	// Free all the data
	free(left);
    free(right);
    free(merged);
}

int main(int argc, char **argv)
{
	int* data;
	int* sorted_data;
	int m, n, id, p, i, output;

	MPI_Status status;

	if (argc != 2) {
		printf("usage: %s <num_items>\n", argv[0]);
		return 1;
	}

	// Number of items to be sorted
	n = atoi(argv[1]);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id); 
	MPI_Comm_size(MPI_COMM_WORLD, &p); // All processors are informed that the # of processors is assigned to p

	if(id==0)
	{
		// data generation
		srandom(0);
		// Make sure that n is a multiple of p
		data = (int *)malloc(n*sizeof(int));
		for(i=0; i<n; i++)
			data[i] = random();
	}


	start_time = clock();

    // TODO: fill in the code here to 
	// (1) distribute the data across the processes
	// (2) sort the data
	// (3) merge sorted data


	if (p == 1) { // When # of processors is 1, no need to use mpi
		mergesort(data, 0, n - 1);
		sorted_data = data;
	}

	else {
		int indiv_data_size = n / p;
		int remainder = n % p;
		int *sendcounts = (int *)malloc(p * sizeof(int)); // Array specifying the number of elements to send
		int *displs = (int *)malloc(p * sizeof(int)); // Array specifying the displacement

		for (i = 0; i < p; i++) {
			sendcounts[i] = indiv_data_size + (i < remainder ? 1 : 0);
			displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
		}

		int local_data_size = sendcounts[id];
		int* local_data = (int *)malloc(local_data_size * sizeof(int));

		// Scatter data to each processes.
		// The reason using Scatterv instead of Scatter is for the case that each processes doesn't have same number of data
		MPI_Scatterv(data, sendcounts, displs, MPI_INT, local_data, local_data_size, MPI_INT, 0, MPI_COMM_WORLD);
		
		// Proceed mergesort for each processes' local data
		mergesort(local_data, 0, local_data_size -1);
 
		int total_participant_num = p;
		int odd_participant_happen = 0; // To check if thre was case that number of participating processes is odd

		// Merge data from processes
		while (total_participant_num > 1) {

			// Total number of participating processes is odd
			if (total_participant_num % 2 != 0) {
				// The middle process send data to 0th process
				if (id == total_participant_num / 2) {
					MPI_Send(local_data, local_data_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
					free(local_data);
					local_data = NULL;
					break;
				}
				// 0th process receive data from the middl process and merge data
				if (id == 0) {
					int opponent_idx = total_participant_num / 2;
					int opponent_data_size = sendcounts[opponent_idx];
					int *opponent_data = (int *)malloc(opponent_data_size * sizeof(int));
					MPI_Recv(opponent_data, opponent_data_size, MPI_INT, opponent_idx, 0, MPI_COMM_WORLD, &status);

					// Merge received data
					int *merged_data = merge(local_data, opponent_data, local_data_size, opponent_data_size);
					free(local_data);
					local_data = merged_data;
					local_data_size += opponent_data_size;

					// Update sendcounts
					sendcounts[id] = local_data_size;
					sendcounts[opponent_idx] = 0;
					free(opponent_data);
				}
				odd_participant_happen = 1; // odd happend -> True
			}
            
			// No matter odd happen or not, we have even number of participating processes
			// First half among processes receive data
			if (id < total_participant_num / 2) { 
				int opponent_idx = total_participant_num - id -1;
				int opponent_data_size = sendcounts[opponent_idx];
				int *opponent_data = (int *)malloc(opponent_data_size * sizeof(int));
				MPI_Recv(opponent_data, opponent_data_size, MPI_INT, opponent_idx, 0, MPI_COMM_WORLD, &status);

				// Merge received data with local data
				int *merged_data = merge(local_data, opponent_data, local_data_size, opponent_data_size);

				free(local_data);

				local_data = merged_data;
				local_data_size += opponent_data_size;

				// Updating sendcounts 
				for (int i = 0; i < total_participant_num / 2; i++) {
					sendcounts[i] += sendcounts[total_participant_num - 1 - i];
					sendcounts[total_participant_num - 1 - i] = 0;
				}
				free(opponent_data);
			}

			else if (id >= total_participant_num / 2){
				// The middle process already sent it's data to 0th process
				if (odd_participant_happen && id == total_participant_num / 2) {
					break;
				}
				else {
					// Send data to opponent
					int opponent_id = total_participant_num - id - 1;
					MPI_Send(local_data, local_data_size, MPI_INT, opponent_id, 0, MPI_COMM_WORLD);
					free(local_data);
					local_data = NULL;
					break;
				}	
			}
				
			// Update step and total_participant_num for next iteration
			total_participant_num = total_participant_num / 2;
			odd_participant_happen = 0;
		} 

	if (id == 0) {sorted_data = local_data;} 
	else {free(local_data);}

	free(sendcounts);
	free(displs);
			          
	}

	stop_time = clock();


	// Check correctness & print execution time
	if (id == 0)
	{
		FILE * fp;

		printf("%d procs, %d items, %f seconds\n", p, n, (stop_time-start_time)/CLOCKS_PER_SEC);

		if(sorted_data==NULL){
			printf("error: sorted_data does not exist\n");
		}
		else{
			for(i = 0; i < n - 1; i++){
				if(sorted_data[i] > sorted_data[i+1])
					printf("error: sorted_data[%d] is greater than sorted_data[%d]\n",i,i+1);
			}
		}
	}
	MPI_Finalize();
}
