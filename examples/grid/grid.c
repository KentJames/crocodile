#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include "grid.h"

// Assume all uvw are zero. This eliminates all coordinate
// calculations and grids only into the middle.
//#define ASSUME_UVW_0

static const double c = 299792458.0;

static inline double uvw_lambda(struct bl_data *bl_data,
                         int time, int freq, int uvw) {
    return bl_data->uvw[3*time+uvw] * bl_data->freq[freq] / c;
}

static inline int coord(int grid_size, double theta,
                 struct bl_data *bl_data,
                 int time, int freq) {
#ifdef ASSUME_UVW_0
    int x = 0, y = 0;
#else
    int x = (int)floor(theta * uvw_lambda(bl_data, time, freq, 0) + .5);
    int y = (int)floor(theta * uvw_lambda(bl_data, time, freq, 1) + .5);
#endif
    return (y+grid_size/2) * grid_size + (x+grid_size/2);
}

static inline void frac_coord(int grid_size, int kernel_size, int oversample,
                       double theta,
                       struct bl_data *bl_data,
                       int time, int freq,
                       int *grid_offset, int *sub_offset) {
#ifdef ASSUME_UVW_0
    double u = 0, v = 0;
#else
    double u = theta * uvw_lambda(bl_data, time, freq, 0);
    double v = theta * uvw_lambda(bl_data, time, freq, 1);
#endif
    int flx = (int)floor(u + .5 / oversample);
    int fly = (int)floor(v + .5 / oversample);
    int xf = (int)floor((u - (double)flx) * oversample + .5);
    int yf = (int)floor((v - (double)fly) * oversample + .5);
    *grid_offset =
        (fly+grid_size/2-kernel_size/2)*grid_size +
        (flx+grid_size/2-kernel_size/2);
    *sub_offset = kernel_size * kernel_size * (yf*oversample + xf);
}



static inline void printkernel(double complex *kernel, int size_x, int size_y){

    int i,j;
    for(j=0;j<size_y;j++){
        for(i=0;i<size_x;i++){
            printf("%.2f+%.2f ",creal(*(kernel+j*size_x + i)),cimag(*(kernel+j*size_x+i)));
        }
        printf("\n");
    }
    printf("\n");
}

static inline void initkernel(double complex *kernel, int size_x, int size_y, double complex init_value){

    int i,j;
    for(j=0;j<size_y;j++){
        for(i=0;i<size_x;i++){
            *(kernel + j*size_x + i) = init_value;
        }
    }


}


void convolve2d_kernels(double complex *kernel1, double complex *kernel2, double complex *output_kernel, int d_in1_x, int d_in1_y, int d_in2_x, int d_in2_y, int d_out_x, int d_out_y){

    //Set middle of output kernel to kernel1. Then we convolve over that kernel with kernel2.
    //
    //
    //printf("Dimensions: %d, %d",d_in1_x,d_in1_y);
    //printf("Still alive! \n");
    int intk_x = d_out_x + (d_in2_x-1);
    int intk_y = d_out_y + (d_in2_y-1);
    double complex *intermediate_kernel = malloc(intk_x * intk_y * sizeof(double complex));
    initkernel(intermediate_kernel,intk_x,intk_y,0+0i);
    int i,j;

    //printf("Still alive! \n");

    for(j=(d_in2_y-1);j<(d_in2_y+(d_in1_y-1));j++){
        for(i=(d_in2_x-1);i<(d_in2_x+(d_in1_x-1));i++){
            *(intermediate_kernel+j*intk_x+i) = *(kernel1+(j-(d_in2_y-1))*d_in1_x+(i-(d_in2_x-1)));
        }
    }
    
    //printf("Intermediate Kernel: \n");
//    printkernel(intermediate_kernel,intk_x,intk_y);

    int ii,jj;
    for(j=0;j<d_out_y;j++){
        for(i=0;i<d_out_x;i++){

            for(jj=0;jj<d_in2_y;jj++){
                for(ii=0;ii<d_in2_x;ii++){
                
                   *(output_kernel+j*d_out_y+i)+= *(intermediate_kernel+(j+jj)*intk_x+(i+ii)) * *(kernel2+jj*d_in2_y+ii);
                }
            }
        }
    }

    free(intermediate_kernel);

}


void weight(unsigned int *wgrid, int grid_size, double theta,
            struct vis_data *vis) {

    // Simple uniform weighting
    int bl, time, freq;
    memset(wgrid, 0, grid_size * grid_size * sizeof(unsigned int));
    for (bl = 0; bl < vis->bl_count; bl++) {
        for (time = 0; time < vis->bl[bl].time_count; time++) {
            for (freq = 0; freq < vis->bl[bl].freq_count; freq++) {
                wgrid[coord(grid_size, theta, &vis->bl[bl], time, freq)]++;
            }
        }
    }
    for (bl = 0; bl < vis->bl_count; bl++) {
        for (time = 0; time < vis->bl[bl].time_count; time++) {
            for (freq = 0; freq < vis->bl[bl].freq_count; freq++) {
                vis->bl[bl].vis[time*vis->bl[bl].freq_count + freq]
                    /= wgrid[coord(grid_size, theta, &vis->bl[bl], time, freq)]; 
                
        }
    }

}}

uint64_t grid_simple(double complex *uvgrid, int grid_size, double theta,
                         struct vis_data *vis) {

    uint64_t flops = 0;
    int bl, time, freq;
    for (bl = 0; bl < vis->bl_count; bl++) {
        for (time = 0; time < vis->bl[bl].time_count; time++) {
            for (freq = 0; freq < vis->bl[bl].freq_count; freq++) {
                uvgrid[coord(grid_size, theta, &vis->bl[bl], time, freq)]
                    += vis->bl[bl].vis[time*vis->bl[bl].freq_count + freq];
                flops += 2;
            }
        }
    }
    return flops;
}

uint64_t grid_wprojection(double complex *uvgrid, int grid_size, double theta,
                          struct vis_data *vis, struct w_kernel_data *wkern) {

    uint64_t flops = 0;
    int bl, time, freq;
#ifdef _ARM_FEATURE_SVE
    //Generate our Predicate Vectors.
    svbool_t pg_t svptrue_b64();
    svbool_t pg_f svpfalse();
    svbool_t pg_ft = svzip1_b64(pg_f, pg_t); 
#endif


#   pragma omp parallel for schedule(guided,100)
    for (bl = 0; bl < vis->bl_count; bl++) {
        for (time = 0; time < vis->bl[bl].time_count; time++) {
            for (freq = 0; freq < vis->bl[bl].freq_count; freq++) {
                // Calculate grid and sub-grid coordinates
                int grid_offset, sub_offset;
                frac_coord(grid_size, wkern->size_x, wkern->oversampling,
                           theta, &vis->bl[bl], time, freq,
                           &grid_offset, &sub_offset);
                // Determine w-kernel to use
                double w = uvw_lambda(&vis->bl[bl], time, freq, 2);
                int w_plane = (int)floor((w - wkern->w_min) / wkern->w_step + .5);
                double complex *wk = wkern->kern_by_w[w_plane].data;
                // Get visibility
                double complex v = vis->bl[bl].vis[time*vis->bl[bl].freq_count+freq];
                // Copy kernel
                int x, y;



#               ifdef _ARM_FEATURE_SVE
                
                uint64_t *igm_r = malloc(sizeof(uint64_t)*vl);
                uint64_t *igm_c = malloc(sizeof(uint64_t)*vl);

                uint64_t *uvgrid_r = malloc(sizeof(uint64_t)*vl);
                uint64_t *uvgrid_c = malloc(sizeof(uint64_t)*vl);

                int index=0;
                
                for (y=0; y < wkern->size_y; y++){

                    for (x = 0; x< wkern->size_x; x++){

                        if(index < (vl-1)){
                            igm_r[index] = sub_offset + y*wkern->size_x + x;
                            igm_c[index] = sub_offset + y*wkern->size_x + x + 1;
                            uvgrid_r[index] = grid_offset + y*wkern->size_x + x;
                            uvgrid_c[index] = grid_offset + y*wkern->size_x + 1;
                            index++;
                        }
                        else { 
                            //If all indexes generated then run the complex multiply!
                            //
                            //Generate indexes in SVE format.
                            svuint64_t ig_r = svld1_u64(pg_t,igm_r);
                            svuint64_t ig_c = svld1_u64(pg_t,igm_c);
                            svuint64_t uv_r = svld1_u64(pg_t,uvgrid_r);
                            svuint64_t uv_c = svld1_u64(pg_t,uvgrid_c);
                            //Load the real and complex parts into their seperate registers.
                            svfloat64_t sub_r = svld1_gather_index(pg_t,(double *)&wk,ig_r);
                            svfloat64_t sub_c = svld1_gather_index(pg_t,(double *)&wk,ig_c);

                            //Now for the actual arithmetic..
                            svfloat64_t r_r_m = svmul_z(pg_t,sub_r,creal(v));
                            svfloat64_t r_c_m = svmul_z(pg_t,sub_r,cimag(v));
                            svfloat64_t c_c_m = svmul_z(pg_t,sub_c,cimag(v));
                            c_c_m = svneg_m(c_c_m,pg_t,c_c_m);
                            svfloat64_t c_r_m = svmul_z(pg_t,sub_c,creal(c));
                            //Add our numbers together
                            r_r_m = svadd_m(pg_t,r_r_m,c_c_m);
                            c_c_m = svadd_m(pg_t,r_c_m,c_r_m);

                            //At long last, add them to the grid.
                            svfloat64_t sub_uv_r = svld1_gather_index(pg_t,(double *)&wk,uv_r);
                            svfloat64_t sub_uv_c = svld1_gather_index(pg_t,(double *)&wk,uv_c);
                            r_r_m = svadd_m(pg_t,sub_uv_r,r_r_m);
                            c_c_m = svadd_m(pg_t,sub_uv_c,c_c_m);
                            svst1_scatter_index(pg_t,(double *)&uvgrid,uv_r,r_r_m);
                            svst1_scatter_index(pg_t,(double *)&uvgrid,uv_c,c_c_m); 

                            index=0;

                        }

                    }

                }

                 

                
#               else
                for (y = 0; y < wkern->size_y; y++) {
                    for (x = 0; x < wkern->size_x; x++) {

                        uvgrid[grid_offset + y*grid_size + x]
                            += v * conj(wk[sub_offset + y*wkern->size_x + x]);
                    }
                }
#               endif
            }
        }
    }
    return flops;
}





void convolve_aw_kernels(struct bl_data *bl,
                         struct w_kernel_data *wkern,
                         struct a_kernel_data *akern) {

    assert(wkern->size_x == akern->size_x);
    assert(wkern->size_y == akern->size_y);
    int size_x = wkern->size_x, size_y = wkern->size_y;
    int ov = wkern->oversampling;

    // We assume that every time channel has their own w-kernel
    const int awkern_count = bl->time_count * akern->freq_count;
    const int awkern_size = size_x * size_y * ov * ov;
    bl->awkern = (double complex *)malloc(awkern_count * awkern_size * sizeof(double complex));

    int time, freq;


//# pragma omp parallel for
    for (time = 0; time < bl->time_count; time++) {
        for (freq = 0; freq < akern->freq_count; freq++) {
            double t = bl->time[time];
            int atime = (int)floor((t - akern->t_min) / akern->t_step + .5);
            int a1i = bl->antenna1 * akern->time_count * akern->freq_count + atime * akern->freq_count + freq;
            int a2i = bl->antenna2 * akern->time_count * akern->freq_count + atime * akern->freq_count + freq;
            struct a_kernel *a1k = &akern->kern_by_atf[a1i];
            struct a_kernel *a2k = &akern->kern_by_atf[a2i];
            double w = bl->uvw[time*3+2] * a1k->freq / c;
            int w_plane = (int)floor((w - wkern->w_min) / wkern->w_step + .5);
            struct w_kernel *wk = &wkern->kern_by_w[w_plane];

            // James Kent: this is my stab at convolving the kernels. Not sure if it is right yet
            // TODO: Implement oversampling.
            int output_size_x = size_x*2 -1;
            int output_size_y = size_y*2 -1;
            double complex *a1a2i = malloc(output_size_x*output_size_y*sizeof(double complex));
            convolve2d_kernels(a1k->data,a2k->data,a1a2i,size_x,size_y,size_x,size_y,output_size_x,output_size_y);

            int output_size_x_awk = output_size_x+(size_x-1);
            int output_size_y_awk = output_size_y+(size_y-1);
            double complex *awk = malloc(output_size_x_awk*output_size_y_awk*sizeof(double complex));
            //printf("Calculating second convolution kernel");
            convolve2d_kernels(a1a2i,wk->data,awk,output_size_x,output_size_y,size_x,size_y,output_size_x_awk,output_size_y_awk);
            //printf("Calculated awkern size: %d",awkern_size);
            //printf("My awkern size!: %d", output_size_x_awk*output_size_y_awk);

            memcpy(&bl->awkern[(time * akern->freq_count + freq) * awkern_size],
                   awk,
                   (output_size_x_awk*output_size_y_awk) * sizeof(double complex));
    
            free(a1a2i);
            free(awk);

        }
    }
}

uint64_t grid_awprojection(double complex *uvgrid, int grid_size, double theta,
                           struct vis_data *vis,
                           struct w_kernel_data *wkern,
                           struct a_kernel_data *akern,
                           int bl_min, int bl_max) {

    // Note that we require awkern to be set on all baselines
    // processed here!

    uint64_t flops = 0;
    int bl, time, freq;
    for (bl = bl_min; bl < bl_max; bl++) {
        const int awkern_size = akern->size_x * akern->size_y *
                                wkern->oversampling * wkern->oversampling;
        const struct bl_data *pbl = &vis->bl[bl];
        for (time = 0; time < pbl->time_count; time++) {
            for (freq = 0; freq < pbl->freq_count; freq++) {
                // Calculate grid and sub-grid coordinates
                int grid_offset, sub_offset;
                frac_coord(grid_size, wkern->size_x, wkern->oversampling,
                           theta, &vis->bl[bl], time, freq,
                           &grid_offset, &sub_offset);
                // Determine kernel frequency, get kernel
                int afreq = (int)floor((pbl->freq[freq] - akern->f_min)
                                       / akern->f_step + .5);
                double complex *awk = &pbl->awkern[
                     (time * akern->freq_count + afreq) * awkern_size];
                // Get visibility
                double complex v = vis->bl[bl].vis[time*vis->bl[bl].freq_count+freq];
                int x, y;
                for (y = 0; y < wkern->size_y; y++) {
                    for (x = 0; x < wkern->size_x; x++) {
                        uvgrid[grid_offset + y*grid_size + x]
                            += v * conj(awk[sub_offset + y*wkern->size_x + x]);
                    }
                }
                flops += 8 * wkern->size_x * wkern->size_y;
            }
        }
    }

    return flops;
}

void make_hermitian(double complex *restrict uvgrid,const int grid_size) {

    struct timeval time1;
    struct timeval time2;

    gettimeofday(&time1,NULL);


//Optimised over reference implementation
#ifdef OPT_OMP
// If running on aarch64 with SVE support, use intrinsics to leverage SVE registers.
#   ifdef _ARM_FEATURE_SVE

        int i_s = 0;
        int lb;
        if (grid_size % 2 == 0) {
            i_s = grid_size + 1;
            lb = ((grid_size*grid_size)-24000)/2; 
        } else {
            i_s = 0;
            lb = (grid_size*grid_size)/2;

        }

        int gs = grid_size*grid_size-1;
        int i; // Loop iterator
        
        svfloat64_t test;
        uint64_t vl = svlen(test);
        printf("Vector Width: %d", vl); 
        svbool_t pg_t = svptrue_b64(); //Used for straight vector load of p0.
        svbool_t pg_f = svpfalse(); //Will interleave this with pg_t to create FNeg predicate.
        svbool_t pg_ft = svzip1_b64(pg_f,pg_t); //Interleaved False/True/False/True.. predicate.

        // Generate indexes for gather/scatter operations.
        uint64_t *igm = malloc(sizeof(uint64_t)*vl);
        for(int i = 0; i < vl - 1; i += 2 ) {
            igm[i] = vl - i - 2;
            igm[i+1] = vl - i - 1;
            printf("Index %d: %d %d \n",i,igm[i],igm[i+1]);

        }
        svuint64_t ig = svld1_u64(pg_t,igm); // Put indexes in SVE format.


        for(i=0;i<(lb-(vl/2)-1);i+=(vl/2)){
            //Load the first sub array of complex numbers. 
            //Cast pointer as a double, as double complex is just interleaved(Re/Im/Re/Im.. etc)  doubles
            svfloat64_t subm1 = svld1(pg_t,(double *)&uvgrid[i]); 
            //Gather Load p1 into vector register using indexes
            svfloat64_t subm2 = svld1_gather_index(pg_t,(double *)&uvgrid[gs-i-(vl/2)],ig); //Gather load.

            //Negate imaginary values and add to other side of grid, then store in memory.
            svfloat64_t subm_neg = svneg_m(subm2, pg_ft,subm2); 
            svfloat64_t subm_add = svadd_m(pg_t,subm1,subm_neg);
            svst1(pg_t, (double*)&uvgrid[i], subm_add);

            //Negate  imaginary values and to other side of grid, then store in memory
            subm_neg = svneg_m(subm1, pg_ft, subm1);
            subm_add = svadd_m(pg_t,subm2,subm_neg);
            svst1_scatter_index(pg_t, (double*)&uvgrid[gs-i-(vl/2)],ig,subm_add); //Scatter store.       
        }
        // Whatever doesn't fit into SVE registers we just iterate through in serial. 
        for(;i<(lb-1);i++){
            double complex g0 = uvgrid[i+i_s];
            uvgrid[i+i_s] = conj(uvgrid[gs-i]);
            uvgrid[gs-i] = conj(g0);

        }
        // Should end up on zero frequency, fingers crossed!
        assert( (i+i_s) == (gs-i) && (i+i_s) == (grid_size+1) * (grid_size/2) );
        uvgrid[(grid_size+1)*(grid_size/2)] += conj(uvgrid[i]);
#   else

        printf("OMP Hermitian Shift \n");
        int i_s = 0;
        int lb;
        if (grid_size % 2 == 0) {
            i_s = grid_size + 1;
            lb = ((grid_size*grid_size)-24000)/2; 
        } else {
            i_s = 0;
            lb = (grid_size*grid_size)/2;

        }

        int gs = grid_size*grid_size-1;
        int i; // Loop iterator
#       pragma omp parallel for
        for(i=0;i<(lb-1);i++){
            double complex g0 = uvgrid[i+i_s];
            uvgrid[i+i_s] += conj(uvgrid[gs-i]);
            uvgrid[gs-i] += conj(g0);
        }

        // Should end up exactly on the zero frequency. Left comments in below for debugging.
    //    printf("is: %d \n",i_s);
    //    printf("Expected end point: %d \n",(grid_size+1)*(grid_size/2));
    //    printf("My end point: %d \n \n",(i+i_s));
    //    printf("gs-i: %d \n \n",gs-i);
    //    assert( (i+i_s) == (gs-i) && (i+i_s) == (grid_size+1) * (grid_size/2) );
        uvgrid[(grid_size+1)*(grid_size/2)] += conj(uvgrid[i]);
#   endif
#else


    printf("Non-OMP Hermitian Shift \n");
    complex double *p0;
    if (grid_size % 2 == 0) {
        p0 = uvgrid + grid_size + 1;
    } else {
        p0 = uvgrid;
    }
    complex double *p1 = uvgrid + grid_size * grid_size - 1;
    // Now simply add cells from the other side of the grid
    while(p0 < p1) {
        double complex g0 = *p0;
        *p0++ += conj(*p1);
        *p1-- += conj(g0);
    }

    // Should end up exactly on the zero frequency
    assert( p0 == p1 && p0 == uvgrid + (grid_size+1) * (grid_size/2) );
    *p0 += conj(*p0);
#endif
    gettimeofday(&time2,NULL);
    printf("Hermitian Time Taken: %f \n", ((double)time2.tv_sec+(double)time2.tv_usec * .000001)-((double)time1.tv_sec+(double)time1.tv_usec * .000001));
}


void fft_shift(double complex *uvgrid, int grid_size) {

    struct timeval time1;
    struct timeval time2;

    gettimeofday(&time1,NULL);
    // Shift the FFT
    assert(grid_size % 2 == 0);
    int x, y;
    for (y = 0; y < grid_size; y++) {
#ifdef OPT_OMP
#   pragma omp parallel for
#endif
        for (x = 0; x < grid_size/2; x++) {
            int ix0 = y * grid_size + x;
            int ix1 = (ix0 + (grid_size+1) * (grid_size/2)) % (grid_size*grid_size);
            double complex temp = uvgrid[ix0];
            uvgrid[ix0] = uvgrid[ix1];
            uvgrid[ix1] = temp;
        }
    }

    gettimeofday(&time2,NULL);
    printf("FFT Shift Time Taken: %f \n", ((double)time2.tv_sec+(double)time2.tv_usec * .000001)-((double)time1.tv_sec+(double)time1.tv_usec * .000001));

}
