#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>
#include <math.h>
#include <complex.h>

#include "grid.h"

double complex predict_dft(double *points, // l/m/l/m/l/m/..
			   double pn,
			   double u,
			   double v,
			   double w){

    double complex vis;
    for(int i = 0; i<pn; ++i){
	double l = *(points + i*2);
	double m = *(points + i*2 + 1);
	double n = sqrt(1.0-l*l-m*m) - 1;
	vis += exp(-2*I * M_PI * (u*l + v*m + w*n));
    }
    return vis;

}
			   
double complex predict_sze(double *points, // l/m/l/m/l/m/..
		       double pn,
		       double u,
		       double v,
		       double w,
		       double du,
		       double dw,
		       int apply_aa,
		       int apply_aa_w,
		       struct sep_kernel_data *kernel_uv,
		       struct sep_kernel_data *kernel_w){

    
    double complex vis;
    for(int i = 0; i<pn; ++i){
	double l = *(points + i*2);
	double m = *(points + i*2 + 1);
	double n = sqrt(1.0-l*l-m*m) - 1;
       	double a = 1.0;
	if (apply_aa) a *= kernel_uv->data[(kernel_uv->size*kernel_uv->oversampling)/2 +
					   (int)(du*kernel_uv->oversampling)];
	if (apply_aa_w) a *= kernel_w->data[(kernel_uv->size*kernel_uv->oversampling)/2 +
	(int)(dw*kernel_uv->oversampling)];
	vis += exp(-2*I * M_PI * (u*l + v*m + w*n));
    }
    return vis;
}

uint64_t predict_grid(double complex *points,
		      double *uvw,
		      double d_u,
		      double d_v,
		      double d_w,
		      struct sep_kernel_data *kernel_uv,
		      struct sep_kernel_data *kernel_w){


}






int main(int argc, char*argv[]){

    int option_index = 0;
    double theta =0, lambda = 0, x0 = 0.25;
    char *sepkernf = NULL, *sepkernwf = NULL;
    int points = 0;
    
    int invalid = 0;
    int dft_compare = 0;


    static struct option options[] =
	{
	 {"theta", required_argument, 0, 't' },
	 {"lambda", required_argument, 0, 'l'},
	 {"sepkern", required_argument, 0, 's'},
	 {"sepkern_w", required_argument, 0, 'w'},
	 {"points", required_argument, 0, 'p'},
	 {"aa_x0", optional_argument, 0, 'x'},
	 {"dft_compare", optional_argument, 0, 'o'},
	 {0, 0, 0, 0}
	};
    int c;
    while ((c = getopt_long(argc, argv, ":", options, &option_index)) != -1) {
        switch(c) {
        case 't': theta = atof(optarg); break;
        case 'l': lambda = atof(optarg); break;
	case 's': sepkernf = optarg; break;
	case 'w': sepkernwf = optarg; break;
	case 'p': points = atoi(optarg); break;
	case 'x': x0 = atof(optarg); break;
	case 'o': dft_compare = 1; break;
        default: invalid = 1; break;
        }
    }
    if (argc == 1) invalid = 1;
    int grid_size = (int)(lambda * theta);
    if (grid_size <= 0){ printf("Invalid grid config!\n"); invalid = 1;}
    
    if (invalid) {
        printf("usage: %s --theta=THETA --lambda=LAM --sepkern=SEPKERN\n", argv[0]);
        printf("\t--points=POINTS [--dft_compare]");
        printf("\n");
        printf("required arguments:\n");
        printf("\t--theta=THETA         Field of view size (in radians)\n");
        printf("\t--lambda=LAM          uv grid size (in wavelenghts)\n");
        printf("\t--sepkern=SEPKERN  seperable kernel u/v file\n");
	printf("\t--sepkern_w=SEPKERNW  seperable kernel w file\n");
        printf("\t--points=IMAGE        number of points to simulate\n");
        printf("optional arguments:\n");
	printf("\t--aa_x0               Sze-Tan Optimised Function Paramter\n");
        printf("\t--dft_compare           compare error vs dft\n");
        return 1;
    }

    

    
    struct sep_kernel_data *sepkern = malloc(sizeof(struct sep_kernel_data));
    struct sep_kernel_data *sepkernw = malloc(sizeof(struct sep_kernel_data));

    printf("Loading Kernel...");
    load_sep_kern(sepkernf,sepkern);
    printf("Loading W Kernel...");
    load_sep_kern(sepkernwf,sepkernw);
    

    



}
