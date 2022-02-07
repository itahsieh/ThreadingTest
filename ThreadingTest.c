#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <omp.h>
#include <pthread.h>
#include <stdlib.h>

const int nThreads = 4;

// the limits of the integral 
double 
x_start = 0.0, 
x_end = M_PI, 
ExactSolution = 2.0;

// number of the integral
size_t N = 1000000000;

double clock_interval(clock_t t_begin, clock_t t_end){
    return (double)((t_end - t_begin)/CLOCKS_PER_SEC);
}

double sine_integral(double lower_limit, double upper_limit, size_t N){
    double dx = (upper_limit - lower_limit) / (double)N;
    double sum = 0.0;
    for (size_t i = 0; i < N; i++){
        double x = lower_limit + ((double)i+0.5) * dx;
        sum += sin(x);
    }
    sum *= dx;
    return sum;
} 

double OMP_sine_integral(double lower_limit, double upper_limit, size_t N, const int nThreads){
    double dx = (upper_limit - lower_limit) / (double)N;
    double sum = 0.0;
    #pragma omp parallel for num_threads(nThreads) reduction(+: sum) 
    for (size_t i = 0; i < N; i++){
        double x = lower_limit + ((double)i+0.5) * dx;
        sum += sin(x);
    }
    sum *= dx;
    return sum;
} 

double *sum_thread;


void *Thread_sine_integral(void *tid_p){
    size_t tid = *((size_t *)tid_p);

    sum_thread[tid] = 0.0;
    double dx_thread = (x_end - x_start)/N;
    
    // integral partitioning
    size_t N_thread = tid < (N % nThreads) ? N/nThreads + 1 : N/nThreads;
    double lower_limit_thread = x_start + (double)tid * (x_end - x_start)/(double)nThreads;

    for (size_t i = 0; i < N_thread; i++){
        double x = lower_limit_thread + ((double)i+0.5) * dx_thread;
        sum_thread[tid] += sin(x);
    }
    /*
    // leap-frog threading
    for (size_t i = 0; i < N; i++){
         if (i % nThreads == tid){
            double x = x_start + ((double)i+0.5) * dx_thread;
            sum_thread[tid] += sin(x);
         }
    }
    */
    return NULL;
}

int main(void){
    clock_t t_start, t_end;
    double sum, error;

    // pre-check
    assert(x_end >= x_start);
#if 0
    // serial job
    t_start = clock();
    sum = sine_integral(x_start, x_end, N);
    t_end = clock();
    
    error = sum - ExactSolution;
    double serial_elapsed = clock_interval(t_start, t_end);
    printf("Serial error:%e Elapsed:%lf s\n", \
    error, serial_elapsed);
#else
    double serial_elapsed = 10.0;
#endif
#if 1
    // OPENMP parallel job
    double omp_start_time = omp_get_wtime();
    sum = OMP_sine_integral(x_start, x_end, N, nThreads);
    double openmp_elapsed = omp_get_wtime()-omp_start_time;
    error = sum - ExactSolution;
    
    printf("OPENMP error:%e Elapsed:%lf ms Efficiency:%f\n", \
    error, openmp_elapsed, serial_elapsed/(openmp_elapsed*nThreads));
#endif
    
    // pthread parallel job
    size_t thread_id[nThreads];
    pthread_t threads[nThreads];
    int sts;
    pthread_attr_t attr;
    
    /* Init thread attribute to joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    sum_thread = malloc(nThreads * sizeof(sum_thread));
    
    time(&t_start);
    /* Create and execute threads */
	for(size_t i = 0; i < nThreads; i++) {
		thread_id[i] = i;
		sts = pthread_create(&threads[i], NULL, Thread_sine_integral, &thread_id[i]);
		assert(sts == 0); /* Temporary measure */
	}
	/* Join threads */
	for(size_t i = 0; i < nThreads; i++) {
		sts = pthread_join(threads[i], NULL);
		assert(sts == 0); /* Temporary measure */
	}
    // sum up the results from every thread
    sum = 0.0;
    for (size_t tid = 0; tid < nThreads; tid++){
        sum += sum_thread[tid];
    }
    sum *= (x_end - x_start) / (double)N;
    time(&t_end);
    double pthread_elapsed = t_end - t_start;
    free(sum_thread);
    error = sum -ExactSolution;
    
    printf("pthread error:%e Elapsed:%lf s Efficiency:%f\n", \
    error, pthread_elapsed, serial_elapsed/(pthread_elapsed*nThreads));

    return 0; 
}
