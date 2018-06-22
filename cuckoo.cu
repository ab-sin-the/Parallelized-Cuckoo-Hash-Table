#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cuda.h>
#include <ctime>

using namespace std;

__device__ int hashing_d (int element, int a, int b, int c, int p, int n){
	return (unsigned int)(a * element + b) % p % n;
}


int hashing (int element, int a, int b, int c, int p, int n){
	return (unsigned int)(a * element + b) % p % n;
}

__global__ void add_two(int p, int n, int N, int t, int* hash_table, int* hash_table2, int* hash_elements, int* func_table, int* a, int* b, int* c){
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	int func1 = 0;
	int func2 = 1;
	int func3 = 2;
	int hash_element = hash_elements[i];
	unsigned int loca1 = hashing_d(hash_element, a[func1], b[func1], c[func1], p, n);
	unsigned int loca2 = hashing_d(hash_element, a[func2], b[func2], c[func2], p, n);
	unsigned int loca3 = hashing_d(hash_element, a[func3], b[func3], c[func3], p, n);
	atomicAdd(&hash_table2[(unsigned int)loca1], 1);
	atomicAdd(&hash_table2[(unsigned int)loca2], 1);
	atomicAdd(&hash_table2[(unsigned int)loca3], 1);
}

__global__ void fix_location(int p, int n, int N, int t, int* hash_table, int* hash_table2, int* hash_elements, int* func_table, int* a, int* b, int* c){
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	int func1 = 0;
	int func2 = 1;
	int func3 = 2;
	int hash_element = hash_elements[i];
	unsigned int loca1 = hashing_d(hash_element, a[func1], b[func1], c[func1], p, n);
	unsigned int loca2 = hashing_d(hash_element, a[func2], b[func2], c[func2], p, n);
	unsigned int loca3 = hashing_d(hash_element, a[func3], b[func3], c[func3], p, n);
	if (hash_table2[(unsigned int)loca1] == 1){
		func_table[i] = 0;
	}
	if (hash_table2[(unsigned int)loca2] == 1){
		func_table[i] = 1;
	}
	if (hash_table2[(unsigned int)loca3] == 1){
		func_table[i] = 2;
	}
}

__global__ void check(int p, int n, int N, int t, int* hash_table, int* hash_elements, int* func_table, int* a, int* b, int* c, int max_count, int *indicator){
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	int curr_func = func_table[i];
	int hash_element = hash_elements[i];
	unsigned int loca = hashing_d(hash_element, a[curr_func], b[curr_func], c[curr_func], p, n);
	if (hash_table[(unsigned int)loca] != hash_element){
		*indicator = -1;
	}
}

__global__ void insert(int p, int n, int N, int t, int* hash_table, int* hash_elements, int* func_table, int* a, int* b, int* c, int max_count){
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	if (i < N){
		int curr_func = func_table[i];
		int hash_element = hash_elements[i];
		unsigned int loca = hashing_d(hash_element, a[curr_func], b[curr_func], c[curr_func], p, n);

		//hash_table[(unsigned int)loca] = hash_element;
		if (hash_table[(unsigned int)loca] != hash_element){
			func_table[i] = (curr_func + 1) % t;
			curr_func = (curr_func + 1) % t;
			loca = hashing_d(hash_element, a[curr_func], b[curr_func], c[curr_func], p, n);
			hash_table[(unsigned int)loca] = hash_element;
		}
	/*for (int j = 0; j < max_count; j++){
		if (hash_table[(unsigned int)loca] != hash_element){
			func_table[i] = (curr_func + 1) % t;
			curr_func = (curr_func + 1) % t;
			loca = ((unsigned int)(a[curr_func] * hash_element + b[curr_func]) % p % n);
			hash_table[(unsigned int)loca] = hash_element;
		}
		__syncthreads();
	}*/

	}

}




__global__ void find(int p, int n, int N, int t, int* hash_table, int* find_elements, int* a, int* b, int* c, int* find_result){
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	int find = 0;
	for (int j = 0; j < t; j++){
		if (hash_table[hashing_d(find_elements[i], a[j], b[j], c[j], p, n)] == find_elements[i]){
			find = 1;
		}
	}
	find_result[i] = find;
}



__global__ void delete_ele(int p, int n, int N, int t, int* hash_table, int* fun_index_table,int* delete_elements, int* a, int* b, int* c, int* find_result){
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	for (int j = 0; j < t; j++){
		if (hash_table[hashing_d(delete_elements[i], a[j], b[j], c[j], p, n)] == delete_elements[i]){
			hash_table[hashing_d(delete_elements[i], a[j], b[j], c[j], p, n)] = 0;
			fun_index_table[hashing_d(delete_elements[i], a[j], b[j], c[j], p, n)] = t;
		}
	}
}

void random_hash_fun(int t, int* a, int* b, int* c){
	for (int i = 0; i < t; i++){
		a[i] = rand();
		b[i] = rand();
		c[i] = rand();
	}
	return;
}

void random_hash_elements(int N, int* hash_elements){
	for (int i = 0; i < N; i++){
		hash_elements[i] = rand() % (1 << 27);
	}
	return;
}

void initialize(int n, int N, int* hash_table, int* hash_table2, int* func_table, int* hash_value, int* hash_elements, int a, int b, int c, int p ){
	for (int i = 0; i < n; i++){
		hash_table[i] = 0;
		hash_table2[i] = 0;
	}
	for (int i = 0; i < N; i++){
		int loca = hashing(hash_elements[i], a, b, c, p, n);
		hash_value[i] = loca;
		func_table[i] = 0;
	}
}



void generating_find_element(int N, int* find_elements, int* hash_element, float partial){
	for (int i = 0; i < (int) (N * partial); i++){
		find_elements[i] = hash_element[i];
	}

	for (int i = (int) (N * partial); i < N; i++){
		find_elements[i] = rand() % (1 << 27);
	}
}

void quicksort(int * hash_value, int* hash_elements, int low, int high)
{
    if(low >= high){
        return;
    }
    int first = low;
    int last = high;
    int key = hash_value[first];
    int key1 = hash_elements[first];
    while(first < last){
        while(first < last && hash_value[last] >= key){
        	--last;
        }

        hash_value[first] = hash_value[last];
        hash_elements[first] = hash_elements[last];
 
        while(first < last && hash_value[first] <= key){
            ++first;
        }
        hash_value[last] = hash_value[first]; 
        hash_elements[last] = hash_elements[first];    
    }
    hash_value[first] = key;
    hash_elements[first] = key1;
    quicksort(hash_value, hash_elements, low, first-1);
    quicksort(hash_value, hash_elements, first+1, high);
}


double once(int p, int n, int t, int N, int max_count, int trail, float partial, int thread_per_block){
	cout << "Trail: " << trail << endl;
	int* hash_table;
	int* hash_table2;
	int* hash_elements;
	int* func_table;
	int* a;
	int* b;
	int* c;
	int* find_result; 
	int* find_elements; 
	int* indicator;
	int count = 0;
	int* hash_value = new int [N];

	cudaMallocManaged(&hash_table, n * sizeof(int));
	cudaMallocManaged(&hash_table2, n * sizeof(int));
	cudaMallocManaged(&hash_elements, N * sizeof(int));
	cudaMallocManaged(&func_table, N * sizeof(int));
	cudaMallocManaged(&a, t * sizeof(int));
	cudaMallocManaged(&b, t * sizeof(int));
	cudaMallocManaged(&c, t * sizeof(int));
	cudaMallocManaged(&find_elements, N * sizeof(int));
	cudaMallocManaged(&find_result, N * sizeof(int));
	cudaMallocManaged(&indicator, sizeof(int));

	*indicator = 0;


	random_hash_fun(t, a, b, c);
	random_hash_elements(N, hash_elements);
	initialize(n, N, hash_table, hash_table2, func_table, hash_value, hash_elements, a[0], b[0], c[0], p);


	quicksort(hash_value, hash_elements, 0, N - 1);

	add_two<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_table2, hash_elements, func_table, a, b, c);
	cudaDeviceSynchronize();
	fix_location<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_table2, hash_elements, func_table, a, b, c);
	cudaDeviceSynchronize();


	for (count = 0; count < max_count; count ++){
		*indicator = 0;
		insert<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_elements, func_table, a, b, c, max_count);
		cudaDeviceSynchronize();
		check<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_elements, func_table, a, b, c, max_count, indicator);
		cudaDeviceSynchronize();
		if (*indicator == 0){
			break;
		}
	}

	check<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_elements, func_table, a, b, c, max_count, indicator);
	cudaDeviceSynchronize();

	int count1 = 0;
	while (*indicator == -1 && count1 < 1000){
		random_hash_fun(t, a, b, c);
		add_two<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_table2, hash_elements, func_table, a, b, c);
		cudaDeviceSynchronize();
		fix_location<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_table2, hash_elements, func_table, a, b, c);
		cudaDeviceSynchronize();

		for (count = 0; count < max_count; count ++){
			*indicator = 0;
			insert<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_elements, func_table, a, b, c, max_count);
			cudaDeviceSynchronize();
			check<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_elements, func_table, a, b, c, max_count, indicator);
			cudaDeviceSynchronize();
			if (*indicator == 0){
				break;
			}
		}
		check<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, hash_elements, func_table, a, b, c, max_count, indicator);
		cudaDeviceSynchronize();
		count1 ++;
	}

	generating_find_element(N, find_elements, hash_elements, partial);

	clock_t start = clock();
	find<<< ceil(N / thread_per_block), thread_per_block>>>(p, n, N, t, hash_table, find_elements, a, b, c, find_result);
	
	cudaDeviceSynchronize();

	double duration = (clock() - start) / (double) CLOCKS_PER_SEC;


	cout << "Time:" << duration << endl;
	cudaMemcpy(hash_table, hash_table, n * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(hash_elements, hash_elements, N * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(find_result, find_result, N * sizeof(int), cudaMemcpyDeviceToHost);

	int val_count = 0;

	val_count = 0;
	for (int i = 0 ; i < N; i++){
		if(find_result[i] == 1){
			val_count ++;
		}
	}

	cout << "Difference between insertion and find result: " << N - val_count << endl;


	cudaFree(hash_table);
	cudaFree(hash_elements);
	cudaFree(find_result);
	cudaFree(a);
	cudaFree(b);
	cudaFree(c);
	cudaFree(hash_table2);
	cudaFree(func_table);
	cudaFree(indicator);
	return duration;
}

int main(void){
	int p = 99984923;
	int count_bit;
	float relative;
	float partial;
	int n;
	int c;
	cout << "Data Bit: " ;
	cin >> count_bit;
	cout << "Find partial: ";
	cin >> partial;
	cout << "Hash table size (0 for 2^25): ";
	cin >> relative; 
	cout << "Evict Chain Constant: ";
	cin >> c;
	int N = (1<<count_bit); 

	if (relative == 0){
		n = 1 << 25;
	}else{
		n = relative * N;
	}

	int t = 3;
	int max_count = (int) c * log(n);
	int thread_per_block;
	if (count_bit >= 18){
		thread_per_block = 1024;
	}else{
		thread_per_block = 64;
	}
	double Time = 0;

	Time += once(p, n, t, N, max_count, 1, partial, thread_per_block);
	Time += once(p, n, t, N, max_count, 2, partial, thread_per_block);
	Time += once(p, n, t, N, max_count, 3, partial, thread_per_block);
	Time += once(p, n, t, N, max_count, 4, partial, thread_per_block);
	Time += once(p, n, t, N, max_count, 5, partial, thread_per_block);

	Time = Time / 5;
	cout << "Average Time: " << Time << endl; 
	cout << "Million of insertion per second: " << (N / Time) / 1000000 << endl;
	return 0;
}
