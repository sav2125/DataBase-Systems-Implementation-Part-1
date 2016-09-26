#include <time.h>
#include<sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include<malloc.h>
#include<limits.h>
#include<math.h>

typedef struct {
	size_t index;
	uint32_t num[625];
} rand32_t;

rand32_t *rand32_init(uint32_t x)
{
	rand32_t *s = malloc(sizeof(rand32_t));
	uint32_t *n = s->num;
	size_t i = 1;
	n[0] = x;
	do {
		x = 0x6c078965 * (x ^ (x >> 30));
		n[i] = x;
	} while (++i != 624);
	s->index = i;
	return s;
}

uint32_t rand32_next(rand32_t *s)
{
	uint32_t x, *n = s->num;
	size_t i = s->index;
	if (i == 624) {
		i = 0;
		do {
			x = (n[i] & 0x80000000) + (n[i + 1] & 0x7fffffff);
			n[i] = (n[i + 397] ^ (x >> 1)) ^ (0x9908b0df & -(x & 1));
		} while (++i != 227);
		n[624] = n[0];
		do {
			x = (n[i] & 0x80000000) + (n[i + 1] & 0x7fffffff);
			n[i] = (n[i - 227] ^ (x >> 1)) ^ (0x9908b0df & -(x & 1));
		} while (++i != 624);
		i = 0;
	}
	x = n[i];
	x ^= (x >> 11);
	x ^= (x <<  7) & 0x9d2c5680;
	x ^= (x << 15) & 0xefc60000;
	x ^= (x >> 18);
	s->index = i + 1;
	return x;
}

int int32_cmp(const void *x, const void *y)
{
	int32_t a = * (const int*) x;
	int32_t b = * (const int*) y;
	return a < b ? -1 : a > b ? 1 : 0;
}

int32_t *generate(size_t n, rand32_t *gen)
{
	size_t i;
	int32_t *a = malloc(n << 2);
	for (i = 0 ; i != n ; ++i)
		a[i] = rand32_next(gen);
	return a;
}

int32_t *generate_sorted_unique(size_t n, rand32_t *gen)
{
	size_t i = 0;
	size_t m = n / 0.7;
	uint8_t z = 0;
	uint32_t *a = malloc(n << 2);
	uint32_t *b = calloc(m, 4);
	while (i != n) {
		uint32_t k = rand32_next(gen);
		if (k != 0) {
			size_t h = (uint32_t) (k * 0x9e3779b1);
			h = (h * (uint64_t) m) >> 32;
			while (b[h] != k) {
				if (b[h] == 0) {
					b[h] = a[i++] = k;
					break;
				}
				if (++h == m) h = 0;
			}
		} else if (z == 0) {
			a[i++] = 0;
			z = 1;
		}
	}
	free(b);
	qsort(a, n, 4, int32_cmp);
	return (int32_t*) a;
}

void ratio_per_bit(const int32_t *a, size_t n)
{
	size_t i, j, *c = calloc(32, sizeof(size_t));
	for (i = 0 ; i != n ; ++i) {
		int32_t x = a[i];
		for (j = 0 ; j != 32 ; ++j)
			c[j] += (a[i] >> j) & 1;
	}
	for (j = 0 ; j != 32 ; ++j)
		fprintf(stderr, "%2ld: %.2f%%\n", j + 1, c[j] * 100.0 / n);
	free(c);
}

int main(int argc, char **argv)
{
	int index2;
	
	rand32_t *gen = rand32_init(time(NULL));
	if (argc < 4)
    	{
       		printf("not enough arguments");
       		exit(0);
    	}
    	int k = atoi(argv[1]);
    	int p = atoi(argv[2]);
	
	int* probes_range = (int*)malloc(sizeof(int)*p);
    	int* number_key_levels = (int*)malloc(sizeof(int)*(argc - 3));
	int* max_key_levels = (int*)malloc(sizeof(int)*(argc - 3));
    	
	
    	int number_key_levels_length = (argc - 3);
    	int min_keys = 0;
    	int max_keys = 0;
    	int multiply_factor = 1;
	int min_multiply_factor = 1;
    	size_t alignment = 32;
	size_t size = 4;

	int32_t * ptr_to_levels[number_key_levels_length];
	int max_keys_single_level = 0;
	int32_t long_max=2147483647;

	//Minimum and maximum number of keys calculated for a given fanout configuration and error printed accordingly
	int j = 3;
    	for(j = 3;j<argc;j++)
    	{
        	number_key_levels[j-3] = atoi(argv[j]);
        	min_keys = j==3 ? 1 : (min_keys + (min_multiply_factor)*(number_key_levels[j-3] - 1));
		min_multiply_factor = j==3 ? 1 : (min_multiply_factor*number_key_levels[j-3]);
        	if(j==3)
            	max_keys_single_level = number_key_levels[j-3] - 1;
        	else
            	max_keys_single_level = multiply_factor*(number_key_levels[j-3] - 1);
		max_key_levels[j-3]=max_keys_single_level;
        	multiply_factor = multiply_factor*number_key_levels[j-3];
        	max_keys = max_keys + max_keys_single_level;
        
		//printf("%d \n",max_keys_single_level);
        	int error = posix_memalign(&ptr_to_levels[j-3], alignment, size*max_keys_single_level);
        	if (error != 0)
           		perror("posix_memalign");
    	}
    	printf("Minimum number of keys required is : %d\n",min_keys);
    	printf("Maximum number of keys required is : %d\n",max_keys);
    	if(k < min_keys || k > max_keys)
    	{
        	printf("error , tree can't be build");
        	exit(0);
    	}
	
	
	//printf("number_key_levels_length %d \n",number_key_levels_length);
	int count=0;
	int32_t *a = generate_sorted_unique(k, gen);
	


	/*int32_t *probes= (int32_t*)malloc(sizeof(int32_t)*(p+1));
	int32_t* a = (int32_t*)malloc(sizeof(int32_t)*k);
	int ppp=0;
	int cnt_p = 0;
	//p=k;
	for(ppp=0;ppp<p-1;ppp++)
	{
		probes[ppp]=cnt_p;
		cnt_p++;
	}
	probes[cnt_p]=8.5;
cnt_p=0;
	for(ppp=0;ppp<k;ppp++)
	{
		a[ppp] = cnt_p ;
		cnt_p++;
	}*/

	

	//printf("Count of max keys in tree :  %d\n",count);
										
	int key_cnt = 0;
	int* count_current_keys = (int*)malloc(sizeof(int)*number_key_levels_length);
	int* flags = (int*)malloc(sizeof(int)*number_key_levels_length);
	int pp = 0 ;
	for(pp=0;pp<number_key_levels_length;pp++)
	{
		flags[pp] = 1;
		count_current_keys[pp] = 0;
		//printf("key is : %d", flags[pp]);
		//printf("count is : %d\n", count_current_keys[pp]);
	}

	int index;
	clock_t begin, end,begin2, end2,begin3, end3;
	double time_spent,time_spent2,time_spent3;

	
	/* here, do your time-consuming job */
	
	begin = clock();
        struct timespec time1,time2;
	double dt2;
	clock_gettime(CLOCK_MONOTONIC,&time1);

	int32_t *probes = generate(p, gen);

	//create tree bottom-up
	while(key_cnt < k)
	{
		
		for(index = number_key_levels_length -1 ; index >= 0; index--)
		{
			if(flags[index] == 1)
			{
				//printf("key_cnt is : %d   ",key_cnt);
				//printf("index is : %d    ",index);
				//printf("count_current_keys[index] : %d\n",count_current_keys[index]);
				ptr_to_levels[index][count_current_keys[index]] = a[key_cnt];
				count_current_keys[index]++;
				int mod_val = count_current_keys[index]%(number_key_levels[index] - 1);
				//printf("mod_val is : %d\n",mod_val);
				if(mod_val == 0)
				{
					flags[index] = 0;
				}
				
				for(index2 = index+1; index2 < number_key_levels_length; index2++)
				{
					flags[index2] = 1;
				}
				break;
			}
			
		}
		key_cnt++;
	}

	//appending int_max to all nodes
	int next_level_nodes=-1;
	for(index = number_key_levels_length -1 ; index >= 0; index--)
	{
		int mod=count_current_keys[index]%(number_key_levels[index] - 1);
		
		//printf("index mod next_level_nodes %d %d %d \n",index,mod,next_level_nodes);
		if(mod!=0)
		{
			int m;
			for(m=0;m<number_key_levels[index]-mod-1;m++)
			{
				ptr_to_levels[index][count_current_keys[index]] = long_max;
				count_current_keys[index]++;
			}
		}	
		if(next_level_nodes!=-1)
		{
			double cur_nodes=(double)((double)next_level_nodes/(double)number_key_levels[index]);
			int current_level_nodes=(int)ceil(cur_nodes);
			if(count_current_keys[index]<current_level_nodes*(number_key_levels[index]-1))
			{
				int m;
				int diff=( current_level_nodes*( number_key_levels[index]-1 )) - count_current_keys[index];
				//printf("1 2 3 %d %d %d\n",(current_level_nodes*(number_key_levels[index]-1)), count_current_keys[index],diff);
				for(m=0 ; m<diff ;m++)
				{
					ptr_to_levels[index][count_current_keys[index]] = long_max;
					count_current_keys[index]++;
				}
			}
		}
		next_level_nodes = count_current_keys[index]/(number_key_levels[index] - 1);
		//printf("next_levek_nodes number_key_levels[index] %d %d \n",next_level_nodes,number_key_levels[index]);

	}

	clock_gettime(CLOCK_MONOTONIC,&time2);
	dt2 = (time2.tv_sec - time1.tv_sec) + (double)(time2.tv_nsec - time1.tv_nsec)*1e-9;
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	
	int32_t key_val = 7;

	begin2 = clock();
	struct timespec t3,t4;
	double dt3;
	clock_gettime(CLOCK_MONOTONIC,&t3);
	//performing binary search to find the range of probes
	for(index2 = 0; index2 < p;index2++)
	{
		int low = 0;
		int high = number_key_levels[0]-2;
		key_val = probes[index2];
		
			//printf("for probe %d: %d %d %d\n",index2,low,high,number_key_levels[0]);
		
		for(index = 0;index < number_key_levels_length;index++)
		{
	
			int result = binary_search(ptr_to_levels[index],low,high,key_val);
			//printf("Binary Search on %d %d \n ",probes[index2],result);
			if(index < number_key_levels_length - 1)
			{
				low = (result) * (number_key_levels[index+1]-1);
				high = low + number_key_levels[index+1]-2;
			}
			else
			{
					int final_result;
					if(result == 0)
						final_result = 0;
					else if(ptr_to_levels[index][result-1] > key_val)
					{
						//printf("probe : %d,arr[result] %d,key_val %d,result : %d\n",index2,ptr_to_levels[index][result-1],key_val,result);
						//printf("inside else if\n");
						final_result = ( (result-1) / ( number_key_levels[index] - 1 ) ) + result - 1;
					}
					else
					{
						//printf("probe : %d,arr[result] %d,key_val %d,result : %d\n",index2,ptr_to_levels[index][result-1],key_val,result);
						//printf("inside else\n");
						final_result = ( (result-1) / ( number_key_levels[index] - 1) ) + result;
					}
					//int final_result =  (result/(number_key_levels[index]-1)) + result + 1;
					probes_range[index2] = final_result;
					//printf("Range Probe %d %d \n",result,index2);
					//printf("forelse probe %d: %d %d %d\n",index2,low,high,number_key_levels[0]);
			}
		}
		//printf("\n");
	}
	clock_gettime(CLOCK_MONOTONIC,&t4);
	dt3 = (t4.tv_sec - t3.tv_sec) + (double)(t4.tv_nsec - t3.tv_nsec)*1e-9;
	end2 = clock();
	
	time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;
	
	//printing tree from leaf level to root level
	printf("Tree: \n");
 	for(index = number_key_levels_length -1 ; index >= 0; index--)
	{
		printf("Level %d\n",index);
		for(index2 = 0 ; index2 < count_current_keys[index]; index2++)
			printf("%d ",ptr_to_levels[index][index2]);
		printf("\n\n");
	}
	
	//Printing ranges corresponding to each probe
	begin3 = clock();
	struct timespec t1,t2;
	clock_gettime(CLOCK_MONOTONIC, &t1);		
	for(index2 = 0; index2 < p;index2++)
	{
		printf("probe is : %d    ",probes[index2]);
		printf("final range is : %d\n",probes_range[index2]);
	}
	clock_gettime(CLOCK_MONOTONIC, &t2);
	double dt = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec)*1e-9;
	end3 = clock();
	time_spent3 = (double)(end3 - begin3) / CLOCKS_PER_SEC;
	
	printf("\n\nTime Phase 1(Creating index and initializing probe): %f seconds  \n", dt2);
	printf("Time Phase 2(Performing probes) : %f seconds \n",dt3);
	printf("Time Phase 3(Printing result to stdout) : %f seconds  \n\n\n",dt);
	
	free(gen);


	int i;
	
	free(number_key_levels);
	free(max_key_levels);
	return EXIT_SUCCESS;
}
	
	
//returns index in a given level passed as array
//low and high specify the index of a specific node
int32_t binary_search(int32_t arr[], int low, int high,int32_t key_val)
{
	//printf("Inside bs low high key_val %d %d %d \n",low,high,key_val);
	
	if(low > high)
	{
		return -1;
	}
	
	/*printf("Inside bs low high key_val %d %d %d \n",low,high,key_val);
	int i=0;	
	for(i=low;i<high;i++)
		printf("%d ",arr[i]);
	printf("\n");*/



	else
	{
		int mid = low + (high - low)/2;
                //printf("%d ffffffffffffff %d fffffff %d\n",mid,arr[low],arr[high]);
		if(key_val < arr[low] && low==0)
			return 0;
		if(key_val < arr[low])
			return low+1;
		else if(key_val >= arr[high])
			return high+1;
		
		if(arr[mid] == key_val)
			return mid+1;
		if(key_val < arr[mid])
		{
			//printf("%d ********* %d\n",mid-1,arr[mid-1]);
			if(mid-1 >= low && key_val >= arr[mid-1])
				return mid;
			binary_search(arr, low, mid-1,key_val);
		}
		else
		{
			if(mid+1 <= high && key_val < arr[mid+1])
				return mid+1;
			binary_search(arr, mid+1, high,key_val);
		}
		
	}

}



