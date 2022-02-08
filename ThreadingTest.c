#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>
#include <omp.h>
#include <pthread.h>
#include <stdlib.h>

// the number of threads
int NThreads,
    NT_init = 1,
    NT_max = 20;

// the limits of the integral 
const double 
XStart = 0.0, 
XEnd = M_PI, 
ExactSolution = 2.0;

// number of the integral
const size_t NInt = 10000000000;

// sum for each pthread's storage
double *sum_thread;
const double dx = (XEnd - XStart)/(double)NInt;
size_t *NInt_thread;

// time interval for clock() function
double clock_interval(clock_t t_begin, clock_t t_end){
    return (double)(t_end - t_begin)/(double)CLOCKS_PER_SEC;
}

// serial sine function integration
double sine_integral(double lower_limit, double upper_limit){
    double sum = 0.0;
    for (size_t i = 0; i < NInt; i++){
        double x = lower_limit + ((double)i+0.5) * dx;
        sum += sin(x);
    }
    sum *= dx;
    return sum;
} 

// OPENMP sine function integration
double sine_integral_OMP(double lower_limit, double upper_limit){
    double sum = 0.0;
    #pragma omp parallel for num_threads(NThreads) reduction(+: sum) 
    for (size_t i = 0; i < NInt; i++){
        double x = lower_limit + ((double)i+0.5) * dx;
        sum += sin(x);
    }
    sum *= dx;
    return sum;
} 

// Threading for sine function integration
void *sine_integral_thread(void *tid_p){
    size_t tid = *((size_t *)tid_p);
    double sum = 0.0;
    double dx_thread = dx;
#if 0
    // integral partitioning by threads  
    double lower_limit_thread[NThreads], NInt_thread[NThreads];
    for (int i = 0; i < NThreads; i++){
        NInt_thread[i] = i < (NInt % NThreads) ? NInt/NThreads + 1 : NInt/NThreads;
        if (i == 0)
            lower_limit_thread[i] = XStart;
        else
            lower_limit_thread[i] = lower_limit_thread[i-1] + (double)NInt_thread[i-1] * dx;
        }
    for (size_t i = 0; i < NInt_thread[tid]; i++){
        double x = lower_limit_thread[tid] + ((double)i + 0.5) * dx_thread;
        sum += sin(x);
    }
#else
    // leap-frog threading
    double XStart_thread = XStart;
    for (size_t i = tid; i < NInt; i += NThreads){
        double x = XStart_thread + ((double)i + 0.5) * dx_thread;
        sum += sin(x);
    }
#endif
    sum_thread[tid] = sum;
    return NULL;
}


/* Main program */
int main(void){
    clock_t t_start, t_end;
    double sum, error, efficiency, speed_up;

    // pre-check: XEnd must be greater than XStart
    assert(XEnd > XStart);
    
    FILE * fptr = fopen("ParallelEfficiency.dat","w");

    #if 1
    // serial job
    fprintf(fptr,"# index 0\n");
    fprintf(fptr,"# Serial: elapsed error\n");
    
    t_start = clock();
    sum = sine_integral(XStart, XEnd);
    t_end = clock();
    
    error = sum - ExactSolution;
    double serial_elapsed = clock_interval(t_start, t_end);
    printf("Serial error:%e Elapsed:%lf s\n", 
        error, serial_elapsed);
    fprintf(fptr,"%lf %e\n", 
        serial_elapsed, error);
    #else
    double serial_elapsed = 10.0;
    #endif

    fprintf(fptr,"\n\n# index 1\n");
    fprintf(fptr,"# OpenMP: NThreads elapsed efficiency error\n");
    for (NThreads = NT_init; NThreads <= NT_max; NThreads++){
        // OpenMP parallel job
        double omp_start_time = omp_get_wtime();
        sum = sine_integral_OMP(XStart, XEnd);
        double openmp_elapsed = omp_get_wtime() - omp_start_time;
        error = sum - ExactSolution;
        efficiency = serial_elapsed/(openmp_elapsed*NThreads);
        speed_up = serial_elapsed/openmp_elapsed;
    
        printf("OPENMP  error:%e Elapsed:%lf s Efficiency:%f\n", 
            error, openmp_elapsed, efficiency);
        fprintf(fptr,"%d %lf %f %e %f\n", 
            NThreads, openmp_elapsed, efficiency, error, speed_up);
    }
    
    fprintf(fptr,"\n\n# index 2\n");
    fprintf(fptr,"# PThread: NThreads elapsed efficiency error\n");
    for(NThreads = NT_init; NThreads <= NT_max; NThreads++){
        // PThread parallel job
        size_t thread_id[NThreads];
        pthread_t threads[NThreads];
        int sts;
        pthread_attr_t attr;
        
        /* Init thread attribute to joinable */
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        
        // PThread routine starts
        time(&t_start);
//         struct timespec ts,te;
//         clock_gettime(CLOCK_MONOTONIC, &ts);
        /* Memory allocation (sum_thread) */
        sum_thread = malloc(NThreads * sizeof(*sum_thread));
        /* Create and execute threads */
        for(size_t i = 0; i < NThreads; i++) {
            thread_id[i] = i;
            sts = pthread_create(&threads[i], NULL, sine_integral_thread, &thread_id[i]);
            /* Temporary measure */
            assert(sts == 0); 
        }
        /* Join threads */
        for(size_t i = 0; i < NThreads; i++) {
            sts = pthread_join(threads[i], NULL);
            /* Temporary measure */
            assert(sts == 0); 
        }
        // sum up the results from every thread
        sum = 0.0;
        for (size_t tid = 0; tid < NThreads; tid++){
            sum += sum_thread[tid];
        }
        sum *= dx;
        time(&t_end);
//         clock_gettime(CLOCK_MONOTONIC, &te);
        free(sum_thread);

        
        double pthread_elapsed = t_end - t_start;
        error = sum -ExactSolution;
        efficiency = serial_elapsed/(pthread_elapsed*NThreads);
        speed_up = serial_elapsed/pthread_elapsed;
        
        printf("PThread error:%e Elapsed:%lf s Efficiency:%f\n", 
            error, pthread_elapsed, efficiency);
        fprintf(fptr,"%d %lf %f %e %f\n", 
            NThreads, pthread_elapsed, efficiency, error, speed_up);
    }
    
    fclose(fptr);
    return 0; 
}
