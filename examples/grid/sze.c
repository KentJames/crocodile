#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>
#include <math.h>
#include <complex.h>
#include <sodium.h>

#include "grid.h"
#include "config.h"


double complex predict_dft(double *points, // l/m/l/m/l/m/..
			   double pn,
			   double u,
			   double v,
			   double w){

    double complex vis = 0 + I*0;
    for(int i = 0; i<pn; ++i){
	double l = *(points + i*2);
	double m = *(points + i*2 + 1);
	double n = sqrt(1.0-l*l-m*m) - 1.0;
	double exponent = -2 * M_PI * (u*l + v*m + w*n);
	vis += cos(exponent) + I*sin(exponent);
    }
    return vis;
}



double complex predict_grid(double *points,
			    int pn,
			    double u,
			    double v,
			    double w,
			    double ov_u,
			    double ov_v,
			    double ov_w,
			    double du,
			    double dw,
			    int aa_support,
			    int aa_support_w,
			    int aa_over,
			    struct sep_kernel_data *grid_conv_uv,
			    struct sep_kernel_data *grid_conv_w,
			    struct sep_kernel_data *grid_corr_lm,
			    struct sep_kernel_data *grid_corr_n){

    // Pre-table our AA kernel.

    double *aa = malloc(pn * sizeof(double));
    for(int p=0; p<pn; ++p){
	double l = points[p*2];
	double m = points[p*2 + 1];
	double n = sqrt(1.0-l*l-m*m) - 1.0;
	
	int lm_size_t = grid_corr_lm->size * grid_corr_lm->oversampling;
	int n_size_t = grid_corr_n->size*grid_corr_n->oversampling;
	double lm_step = 1.0/(double)lm_size_t;
	double n_step = 1.0/(double)n_size_t; // Not too sure about this

	double a = 1.0;
	a *= grid_corr_lm->data[(int)round(((du*l)/lm_step) + (lm_size_t)/2)];
	a *= grid_corr_lm->data[(int)round(((du*m)/lm_step) + (lm_size_t)/2)];
	a *= grid_corr_n->data[(int)round(((dw*n)/n_step) + n_size_t/2)];

	aa[p] = a;
	printf("AA: %f \n",a);
    }


    double complex visg = 0 + I*0; 
    for(int ius = 0; ius < aa_support; ++ius){
	for(int ivs = 0; ivs < aa_support; ++ivs){
	    for(int iws = 0; iws < aa_support_w; ++iws){
		double dus = u + du * ((double)ius - (floor((double)aa_support/2.0)) + ov_u/(double)aa_over);
		double dvs = v + du * ((double)ivs - (floor((double)aa_support/2.0)) + ov_v/(double)aa_over);
		double dws = w + dw * ((double)iws - (floor((double)aa_support_w/2.0)) + ov_w/(double)aa_over);
		printf("dus: %f",dus);
		int aas_u = ius * aa_over + ov_u;
		int aas_v = ivs * aa_over + ov_v;
		int aas_w = iws * aa_over + ov_w;
		double gcf = grid_conv_uv->data[aas_u] * grid_conv_uv->data[aas_v] * grid_conv_w->data[aas_w];
		
		double complex vis = 0 + I*0;
		for(int p=0; p<pn; ++p){
		    double l = *(points + p*2);
		    double m = *(points + p*2 + 1);
		    double n = sqrt(1.0-l*l-m*m) - 1.0;

		    // Now calculate the DFT

		    double exponent = -2 * M_PI * (dus*l + dvs*m + dws*n);
		    double complex visl = cos(exponent) + I * sin(exponent);
		    vis += visl/aa[p];    
		}

		
		double complex visgcf = vis * gcf;
		printf("GCF: %.15f\n",gcf);
		//q		printf("GCF: %.15f + i%.15f\n",creal(visgcf),cimag(visgcf));
		
		visg += visgcf;
	    }
	}
    }
    free(aa);
    return visg;
}


void generate_random_points(double *points,
			    int pn,
			    double theta){

    for(int i=0; i<pn; ++i){
	double l = (double)(randombytes_uniform(65536))/65536 - 0.5;
	double m = (double)(randombytes_uniform(65536))/65536 - 0.5;
	l *= theta;
	m *= theta;
	points[i*2] = l;
	points[i*2 + 1] = m;
    }
}




int main(int argc, char*argv[]){

    int option_index = 0;
    double theta =0, lambda = 0, x0 = 0.25;
    char *sepkernf = NULL, *sepkernwf = NULL,
	*sepkern_imf = NULL, *sepkern_im_wf = NULL;
    int points = 0;
    
    int invalid = 0;
    int dft_compare = 0;


    static struct option options[] =
	{
	 {"theta", required_argument, 0, 't' },
	 {"lambda", required_argument, 0, 'l'},
	 {"sepkern", required_argument, 0, 's'},
	 {"sepkern_w", required_argument, 0, 'w'},
	 {"sepkern_im", required_argument, 0, 'i'},
	 {"sepkern_im_w", required_argument, 0, 'm'},
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
	case 'i': sepkern_imf = optarg; break;
	case 'm': sepkern_im_wf = optarg; break;
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
        printf("\t--theta=THETA                 Field of view size (in radians)\n");
        printf("\t--lambda=LAM                  uv grid size (in wavelenghts)\n");
        printf("\t--sepkern=SEPKERN             seperable kernel u/v file\n");
	printf("\t--sepkern_w=SEPKERNW          seperable kernel w file\n");
	printf("\t--sepkern_im=SEPKERN_IM       seperable kernel l/m file\n");
	printf("\t--sepkern_im_w=SEPKERN_IM_W   seperable kernel n file \n");
        printf("\t--points=IMAGE        number of points to simulate\n");
        printf("optional arguments:\n");
	printf("\t--aa_x0               Sze-Tan Optimised Function Paramter\n");
        printf("\t--dft_compare           compare error vs dft\n");
        return 1;
    }

    struct sep_kernel_data *sepkern = malloc(sizeof(struct sep_kernel_data));
    struct sep_kernel_data *sepkernw = malloc(sizeof(struct sep_kernel_data));
    struct sep_kernel_data *sepkern_im = malloc(sizeof(struct sep_kernel_data));
    struct sep_kernel_data *sepkern_im_w = malloc(sizeof(struct sep_kernel_data));
    
    
    printf("Loading Kernel...");
    load_sep_kern(sepkernf,sepkern);
    printf("Loading W Kernel...");
    load_sep_kern(sepkernwf,sepkernw);

    printf("Loading AA Kernel...");
    load_sep_kern(sepkern_imf, sepkern_im);
    printf("Loading AA Kernel...");
    load_sep_kern(sepkern_im_wf, sepkern_im_w);

    printf("Generating Random Points...");
    double *randompoints = malloc(2 * points*sizeof(double));
    generate_random_points(randompoints,points,theta);    
    printf("done\n");
    
    double complex vis_sze = predict_grid(randompoints,
					  points,
					  1.0,
					  1.0,
					  1.0,
					  43.0, 43.0, 43.0,
					  5.0, 49.9374216791,
					  sepkern->size,
					  sepkernw->size,
					  4096,
					  sepkern,
					  sepkernw,
					  sepkern_im,
					  sepkern_im_w);
    double complex vis = predict_dft(randompoints, points, 1.0, 1.0, 1.0);
    printf("\nVis: Re: %f, Im: %f \n", creal(vis), cimag(vis));				      
    printf("Vis Sze: Re: %f, Im: %f \n", creal(vis_sze), cimag(vis_sze));
    printf("Delta: Re: %f, Im: %f \n", fabs(creal(vis)-creal(vis_sze)),fabs(cimag(vis)-cimag(vis_sze)));
}
