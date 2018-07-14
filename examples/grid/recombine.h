
#ifndef RECOMBINE_H
#define RECOMBINE_H

#include <stdbool.h>
#include <stdint.h>
#include <complex.h>
#include <fftw3.h>

void *read_dump(int size, char *name, ...);
int write_dump(void *data, int size, char *name, ...);

// Generic 1D kernels
double *generate_Fb(int yN_size, int yB_size, double *pswf);
double *generate_Fn(int yN_size, int xM_yN_size, double *pswf);
double *generate_m(int image_size, int yP_size, int yN_size, int xM_size, int xMxN_yP_size,
                   double *pswf);
void prepare_facet(int yB_size, int yP_size,
                   double *Fb,
                   double complex *facet, int facet_stride,
                   double complex *BF, int BF_stride);
void extract_subgrid(int yP_size, int xM_yP_size, int xMxN_yP_size, int xM_yN_size, int subgrid_offset,
                     double *m_trunc, double *Fn,
                     complex double *BF, int BF_stride,
                     complex double *MBF, fftw_plan MBF_plan,
                     complex double *NMBF, int NMBF_stride);
void add_facet(int xM_size, int xM_yN_size, int facet_offset,
               complex double *NMBF, int NMBF_stride,
               complex double *out, int out_stride);

// Configuration for 2D recombination
struct recombine2d_config {
    // Local configuration (might depend on rank)
    char *facet_file; // Facet (debugging)
    char *stream_check; // Check stream contents (debugging)
    double stream_check_threshold;
    char *stream_dump; // File dump I/O (debugging)
    // Sizings
    int image_size; // total image size
    int yB_size; // facet size
    int yN_size; //  ", padded to avoid inaccurate regions of PSWF
    int yP_size; //  ", further padded to minimal size of accurate "m" mask multiplication
    int xA_size; // sub-grid size
    int xM_size; //  ", padded to prevent PSWF convolution from aliasing
    int xMxN_yP_size; // subgrid/facet size, with extra padding for "m" to allow us to work within yP limit
    // Derived sizes
    int xA_yP_size; // Distance between facets (TODO: replace by more flexible facets)
    int xM_yP_size; // Buffer size for "m" mask multiplication
    int xM_yN_size; // Size of subgrid/facet pieces to exchange
    // Derived data layout (F, b*F, m(b*F), n*m(b*F) respectively for two axes)
    size_t F_size, BF_size, MBF_size;
    size_t NMBF_size, NMBF_BF_size, NMBF_NMBF_size;
    size_t SG_size; // subgrid size
    uint64_t F_stride0, F_stride1;
    uint64_t BF_stride0, BF_stride1;
    uint64_t NMBF_stride0, NMBF_stride1;
    uint64_t NMBF_BF_stride0, NMBF_BF_stride1;
    uint64_t NMBF_NMBF_stride0, NMBF_NMBF_stride1;
    // Working arrays (F[b], F[n], m)
    double *Fb, *Fn, *m;
};

bool recombine2d_set_config(struct recombine2d_config *cfg,
                            int image_size, char *pswf_file,
                            int yB_size, int yN_size, int yP_size,
                            int xA_size, int xM_size, int xMxN_yP_size);
void recombine2d_free(struct recombine2d_config *cfg);

uint64_t recombine2d_global_memory(struct recombine2d_config *cfg);
uint64_t recombine2d_worker_memory(struct recombine2d_config *cfg);

struct recombine2d_worker {

    // Configuration
    struct recombine2d_config *cfg;

    // Time (in s) spent in different stages
    double pf1_time, es1_time, ft1_time;
    double pf2_time, es2_time, ft2_time;

    // Plans associated with buffers
    int BF_batch; fftw_plan BF_plan; // shared
    fftw_plan NMBF_BF_plan, MBF_plan;

    // Private buffers
    double complex *MBF;
    double complex *NMBF;
    double complex *NMBF_BF;

};

fftw_plan recombine2d_bf_plan(struct recombine2d_config *cfg, int BF_batch, double complex *BF);
void recombine2d_init_worker(struct recombine2d_worker *worker, struct recombine2d_config *cfg,
                             int BF_batch, fftw_plan BF_plan);
void recombine2d_free_worker(struct recombine2d_worker *worker);

void recombine2d_pf0_ft0_omp(struct recombine2d_worker *worker,
                             complex double *F, complex double *BF);
void recombine2d_es0_pf1_ft1(struct recombine2d_worker *worker,
                             int i0, complex double *BF);
void recombine2d_es1(struct recombine2d_worker *worker,
                     int i0, int i1, double complex *NMBF_NMBF);

void recombine2d_af0_af1(struct recombine2d_config *cfg,
                         double complex *subgrid,
                         int facet_off0, int facet_off1,
                         double complex *NMBF_NMBF);

#endif // RECOMBINE_H