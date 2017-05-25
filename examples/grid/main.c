
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <stdint.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <complex.h>
#include <fftw3.h>
#include <omp.h>

#include "grid.h"

//Returns the baseline distance. Used for qsort. 
int uv_delta(const void * a, const void * b){

    struct bl_data *baselineA = (struct bl_data *)a;
    struct bl_data *baselineB = (struct bl_data *)b;

    return( (sqrt(pow(baselineB->uvw[0],2) + pow(baselineB->uvw[1],2))) 
            - (sqrt(pow(baselineA->uvw[0],2) + pow(baselineA->uvw[1],2)))   ); //sorts longest baselines first
    //return( (sqrt(pow(baselineA->uvw[0],2) + pow(baselineA->uvw[1],2))) 
    //- (sqrt(pow(baselineB->uvw[0],2) + pow(baselineB->uvw[1],2)))   ); //qsorts shortest baselines first
}


int main(int argc, char *argv[]) {
    
    struct timeval time1;
    struct timeval time2;

    gettimeofday(&time1,NULL);

    //int threads = omp_get_num_procs();


    // Read parameters
    static struct option options[] =
      {
        {"theta",   required_argument, 0, 't' },
        {"lambda",  required_argument, 0, 'l' },
        {"wkern",   optional_argument, 0, 'w' },
        {"akern",   optional_argument, 0, 'a' },
        {"grid",    optional_argument, 0, 'g' },
        {"image",   optional_argument, 0, 'i' },
        {"min-bl",  optional_argument, 0, 'b' },
        {"max-bl",  optional_argument, 0, 'B' },
        {"threads", optional_argument, 0, 'P' },
        {0, 0, 0, 0}
      };
    int option_index = 0;
    double theta = 0, lambda = 0;
    char *wkern_file = NULL, *akern_file = NULL,
         *grid_file = NULL, *image_file = NULL;
    double bl_min = DBL_MIN, bl_max = DBL_MAX;
    int threads=1;
    int c; int invalid = 0;
    while ((c = getopt_long(argc, argv, ":", options, &option_index)) != -1) {
        switch(c) {
        case 't': theta = atof(optarg); break;
        case 'l': lambda = atof(optarg); break;
        case 'w': wkern_file = optarg; break;
        case 'a': akern_file = optarg; break;
        case 'g': grid_file = optarg; break;
        case 'i': image_file = optarg; break;
        case 'b': bl_min = atof(optarg); break;
        case 'B': bl_max = atof(optarg); break;
        case 'P': threads = atol(optarg); break;
        default: invalid = 1; break;
        }
    }

    // Check grid parameters
    int grid_size = (int)(theta * lambda);
    size_t grid_byte_size = grid_size * grid_size * sizeof(double complex);
    if (grid_size <= 0) {
        fprintf(stderr, "Invalid grid configuration!\n");
        invalid = 1;
    }

    // Must have an input file
    const char *vis_file = 0;
    if (optind + 1 == argc) {
        vis_file = argv[optind];
    } else {
        printf("Please supply a visibility input file!\n");
        invalid = 1;
    }
    if (invalid) {
        printf("usage: %s --theta=THETA --lambda=LAM [--grid=GRID]\n", argv[0]);
        printf("              [--image=IMAGE] [--wkern=WKERN] [--akern=AKERN]\n");
        printf("              [--min-bl=MIN_BL] [--max-bl=MAX_BL]\n");
        printf("              INPUT\n");
        printf("\n");
        printf("optional arguments:\n");
        printf("  --theta=THETA         Field of view size (in radians)\n");
        printf("  --lambda=LAM          uv grid size (in wavelenghts)\n");
        printf("  --grid=GRID           grid output file\n");
        printf("  --image=IMAGE         image output file\n");
        printf("  --wkern=WKERN         w-kernel file to use for w-projection\n");
        printf("  --akern=AKERN         A-kernel file to use for w-projection\n");
        printf("  --min-bl=MIN_BL       Minimum baseline length to consider (in km)\n");
        printf("  --max-bl=MAX_BL       Maximum baseline length to consider (in km)\n");
        printf("  --threads=THREADS     Number of threads to parallelise for\n");
        printf("positional arguments:\n");
        printf("  input                 input visibilities\n");
        return 1;
    }

    // Intialise HDF5
    init_dtype_cpx();

    // Open files
    struct vis_data vis;
    struct w_kernel_data wkern;
    struct a_kernel_data akern;
    int grid_fd = -1, image_fd = -1;
    if (load_vis(vis_file, &vis, bl_min, bl_max)) {
        return 1;
    }
    if (wkern_file) {
        if (load_wkern(wkern_file, theta, &wkern)) {
            return 1;
        }
    }
    if (wkern_file && akern_file) {
        if (load_akern(akern_file, theta, &akern)) {
            return 1;
        }
    }
    if (grid_file) {
        grid_fd = open(grid_file, O_CREAT | O_TRUNC | O_WRONLY, 0777);
        if (grid_fd == -1) {
            perror("Failed to open grid file");
            return 1;
        }
    }
    if (image_file) {
        image_fd = open(image_file, O_CREAT | O_TRUNC | O_WRONLY, 0777);
        if (image_fd == -1) {
            perror("Failed to open image file");
            return 1;
        }
    }

    // Allocate grid
    printf("\nGrid size:    %d x %d (%.2f GB)\n", grid_size, grid_size, (double)(grid_byte_size)/1000000000);
    double complex *uvgrid = (double complex *)calloc(grid_byte_size, 1); // Our parent grid.
    double complex **uvgrid_cp=calloc(threads,sizeof(double complex **)); //Points to our individual arrays.

    int i;
    for(i=0;i<threads;i++){
        uvgrid_cp[i] = calloc(grid_byte_size,1); //Allocate each subgrid.
    }
    // Sort visibilities in terms of baseline length, which is directly proportional to compute cost.
    printf("Sorting visibilities... \n");
    qsort(vis.bl,vis.bl_count,sizeof(struct bl_data),uv_delta);
    
    // Simple uniform weight (we re-use the grid to save an allocation)
    // We weight with a single thread..
    printf("Weighting...\n");
    weight((unsigned int *)uvgrid, grid_size, theta, &vis);
    memset(uvgrid, 0, grid_size * grid_size * sizeof(unsigned int));

    int bl_ind[threads]; //Indices for dividing the visibility array
    int bl_per_core_t[threads]; //BL per core


   /* // Load spreading attempt one: Using 2Vn/(threads+n) : Not good enough. 
    for(i=0;i<threads;i++){
        bl_ind[i] = (2 * vis.bl_count * i) / (threads + i);
    }
    for(i=0;i<threads;i++){
        if(i<threads-1){
            bl_per_core_t[i] = bl_ind[i+1] - bl_ind[i];
        }  
        else {
            bl_per_core_t[i] = vis.bl_count - bl_ind[i] - 1;
        }
    }

    */


    
    // Load spreading attempt two: Using the fibonacci sequence. Works quite well on 4-8 threads.
    if(threads==1){
        bl_per_core_t[0] = vis.bl_count;
    }
    else{

        int sum;
        int fibonacci[threads];
        fibonacci[0] = 1;
        fibonacci[1] = 3;
        sum = fibonacci[0]+fibonacci[1];
        for(i=2;i<threads;i++){
            fibonacci[i] = fibonacci[i-1] + fibonacci[i-2];
            sum += fibonacci[i];
        }

        int bl_agg = vis.bl_count / sum;
        int bl_rem = vis.bl_count % sum;
        //printf("BL Agg: %d Sum Fib: %d\n",bl_agg,sum);


        for(i=0;i<threads;i++){   
            if(i<threads-1){
                bl_per_core_t[i] = bl_agg * fibonacci[i];
            } 
            else{
                bl_per_core_t[i] = bl_agg * fibonacci[i] + bl_rem;
            }
            //printf("BL Per Core: %d\n",bl_per_core_t[i]);
        }
        int agg = 0;
        for(i=0;i<threads;i++){
            bl_ind[i] = agg;
            agg += bl_per_core_t[i];
            //printf("Indice: %d\n",bl_ind[i]);
        }
    }


    //Lets split up our visibilities to distribute amongst the threads..
    struct vis_data *vis_cp=calloc(threads,sizeof(struct vis_data));

    for(i=0;i<(threads);i++){
        vis_cp[i].bl = calloc(bl_per_core_t[i],sizeof(struct bl_data));
        vis_cp[i].bl_count = bl_per_core_t[i];
        memcpy(vis_cp[i].bl,&vis.bl[bl_ind[i]],bl_per_core_t[i] * sizeof(struct bl_data));
    }


    struct timeval proj1,proj2,proj3;

    // Set up performance counters
    struct perf_counters counters;
    open_perf_counters(&counters);

    uint64_t flops = 0, mem = 0;
    if (!wkern_file) {
        printf("Gridder: Simple imaging\n");
        enable_perf_counters(&counters);
        flops = grid_simple(uvgrid, grid_size, theta, &vis);
        disable_perf_counters(&counters);
        // Assuming 0.5 flop/B
        mem = flops * 2;
    } else if (!akern_file) {

        gettimeofday(&proj1,NULL);
        printf("Gridder: W-projection\n");
        enable_perf_counters(&counters);
#       pragma omp parallel for
        for(i=0;i<threads;i++){
            flops += grid_wprojection(uvgrid_cp[i], grid_size, theta, &vis_cp[i], &wkern);
            printf("Thread %d completed.\n",i);
        }
        gettimeofday(&proj2,NULL);
#       pragma omp parallel for
        for(i=0;i<threads;i++){
            int j;
            for(j=0;j<(grid_size*grid_size);j++){
                uvgrid[j] += *(*(uvgrid_cp+i)+j);
            }
        }

        for(i=0;i<threads;i++){
            free(uvgrid_cp[i]);
        }
        free(vis_cp);
        free(uvgrid_cp);
        
        disable_perf_counters(&counters);
        gettimeofday(&proj3,NULL);

        printf("\n Projection Time taken: %f \n",((double)proj2.tv_sec+(double)proj2.tv_usec * .000001)-((double)proj1.tv_sec+(double)proj1.tv_usec * .000001));
        printf("\n Reduction Time taken: %f \n",((double)proj3.tv_sec+(double)proj3.tv_usec * .000001)-((double)proj2.tv_sec+(double)proj2.tv_usec * .000001));
        printf("\n Total time taken: %f \n",((double)proj3.tv_sec+(double)proj3.tv_usec * .000001)-((double)proj1.tv_sec+(double)proj1.tv_usec * .000001));
        // Assuming 10 flop/B
        mem = flops / 10;
    } else {
        printf("Gridder: AW-projection (ignoring convolution)\n");

        // We need to chunk the convolution so we don't run out of
        // memory...
        const int bl_chunk = 1000;
        int bl_min;
        for (bl_min = 0; bl_min < vis.bl_count; bl_min+=bl_chunk) {
            int bl_max = bl_min + bl_chunk;
            if (bl_max > vis.bl_count) bl_max = vis.bl_count;

            printf("\"Convolving\" %d-%d...\n", bl_min, bl_max-1);
            int bl;
            for (bl = bl_min; bl < bl_max; bl++) {
                convolve_aw_kernels(&vis.bl[bl], &wkern, &akern);
            }            

            // Do convolution
            enable_perf_counters(&counters);
            
            flops += grid_awprojection(uvgrid, grid_size, theta, &vis, &wkern, &akern, bl_min, bl_max);
            disable_perf_counters(&counters);

            // Free convolved kernels
            for (bl = bl_min; bl < bl_max; bl++) {
                free(vis.bl[bl].awkern);
                vis.bl[bl].awkern = NULL;
            }
        }

        // Assuming 0.5 flop/B
        mem = flops * 2;
    }

    // Show performance counters after gridding
    printf("\nCounters:\n");
    print_perf_counters(&counters, flops, mem);

    // Make hermitian
    printf("\nMake hermitian...\n");
    make_hermitian(uvgrid, grid_size);

    // Write grid
    if (grid_fd != -1) {
        printf("Write grid...\n");
        int i;
        for (i = 0; i < grid_size; i++) {
            write(grid_fd, uvgrid + i * grid_size, grid_byte_size / grid_size);
        }
        close(grid_fd);
    }
    if (image_fd != -1) {
        printf("FFT...\n");

        // First shift zero frequency
        fft_shift(uvgrid, grid_size);

        // Do DFT. Complex-to-complex to keep with numpy (TODO: optimize)
        // TODO: Use parallel fftw call. 
        //
        fftw_init_threads();
        fftw_plan plan;
        fftw_plan_with_nthreads(threads);
        plan = fftw_plan_dft_2d(grid_size, grid_size, uvgrid, uvgrid, -1, FFTW_ESTIMATE);
        fftw_execute_dft(plan, uvgrid, uvgrid);

        // Shift zero frequency back into centre
        fft_shift(uvgrid, grid_size);

        // Write real part to disk
        printf("Write image...\n");
        int i;
        double *row = malloc(sizeof(double) * grid_size);
        for (i = 0; i < grid_size; i++) {
            int j;
            for (j = 0; j < grid_size; j++) {
                row[j] = creal(uvgrid[i*grid_size+j]);
            }
            write(image_fd, row, sizeof(double) * grid_size);
        }
        close(image_fd);
    }
    gettimeofday(&time2,NULL);

    printf("\n Time taken: %f \n",((double)time2.tv_sec+(double)time2.tv_usec * .000001)-((double)time1.tv_sec+(double)time1.tv_usec * .000001));
    return 0;
}
