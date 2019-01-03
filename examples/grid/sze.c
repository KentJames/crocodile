#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>
#include <math.h>
#include <complex.h>
#include <sodium.h>
#include "grid.h"

double grid_correction[] = { 0.04951909,  0.04635579,  0.04712115,  0.05060831,  0.0560945,   0.06315665,
			     0.0715563,   0.08116772,  0.09193298,  0.10383379,  0.11687376,  0.13106728,
			     0.14643229,  0.16298567,  0.18074009,  0.19970191,  0.21986963,  0.24123281,
			     0.26377116,  0.28745395,  0.31223949,  0.33807485,  0.36489561,  0.39262584,
			     0.42117816,  0.45045398,  0.48034386,  0.51072799,  0.54147689,  0.57245214,
			     0.6035073,   0.634489,    0.66523804,  0.69559069,  0.72538003,  0.75443738,
			     0.78259381,  0.80968166,  0.83553611,  0.85999675,  0.88290916,  0.90412638,
			     0.92351044,  0.94093371,  0.95628025,  0.96944698,  0.98034476,  0.98889936,
			     0.99505217,  0.9987609,   1.,          0.9987609,   0.99505217,  0.98889936,
			     0.98034476,  0.96944698,  0.95628025,  0.94093371,  0.92351044,  0.90412638,
			     0.88290916,  0.85999675,  0.83553611,  0.80968166,  0.78259381,  0.75443738,
			     0.72538003,  0.69559069,  0.66523804,  0.634489,    0.6035073,   0.57245214,
			     0.54147689,  0.51072799,  0.48034386,  0.45045398,  0.42117816,  0.39262584,
			     0.36489561,  0.33807485,  0.31223949,  0.28745395,  0.26377116,  0.24123281,
			     0.21986963,  0.19970191,  0.18074009,  0.16298567,  0.14643229,  0.13106728,
			     0.11687376,  0.10383379,  0.09193298,  0.08116772,  0.0715563,   0.06315665,
			     0.0560945,   0.05060831,  0.04712115,  0.04635579,  0.04951909};

double grid_correction_w[] = { 0.22637738,  0.2126835,   0.20596314,  0.20474932,  0.2079147,   0.21459257,
			       0.22411589,  0.23597038,  0.24975825,  0.26517025,  0.28196406,  0.29994765,
			       0.3189664,   0.33889324,  0.35962107,  0.38105692,  0.40311757,  0.42572616,
			       0.44880976,  0.47229746,  0.49611909,  0.52020422,  0.54448159,  0.56887872,
			       0.59332168,  0.61773508,  0.64204216,  0.6661649,   0.69002432,  0.71354073,
			       0.73663407,  0.75922431,  0.78123181,  0.80257772,  0.82318442,  0.84297595,
			       0.86187838,  0.8798203,   0.89673317,  0.91255173,  0.92721442,  0.94066368,
			       0.95284636,  0.96371397,  0.97322303,  0.98133529,  0.98801799,  0.99324406,
			       0.99699224,  0.99924729,  0.99999999,  0.99924729,  0.99699224,  0.99324406,
			       0.98801799,  0.98133529,  0.97322303,  0.96371397,  0.95284636,  0.94066368,
			       0.92721442,  0.91255173,  0.89673317,  0.8798203,   0.86187838,  0.84297595,
			       0.82318442,  0.80257772,  0.78123181,  0.75922431,  0.73663407,  0.71354073,
			       0.69002432,  0.6661649,   0.64204216,  0.61773508,  0.59332168,  0.56887872,
			       0.54448159,  0.52020422,  0.49611909,  0.47229746,  0.44880976,  0.42572616,
			       0.40311757,  0.38105692,  0.35962107,  0.33889324,  0.3189664,   0.29994765,
			       0.28196406,  0.26517025,  0.24975825,  0.23597038,  0.22411589,  0.21459257,
			       0.2079147,   0.20474932,  0.20596314,  0.2126835,   0.22637738};

double complex predict_dft(double *points, // l/m/l/m/l/m/..
			   double pn,
			   double u,
			   double v,
			   double w){

    double complex vis = 0 + I*0;
    for(int i = 0; i<pn; ++i){
	double l = *(points + i*2);
	double m = *(points + i*2 + 1);
	double n = sqrt(1.0-l*l-m*m) - 1;
	double exponent = -2 * M_PI * (u*l + v*m + w*n);
	vis += cos(exponent) + I*sin(exponent);
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
			   double theta,
			   int apply_aa,
			   int apply_aa_w,
			   struct sep_kernel_data *kernel_uv,
			   struct sep_kernel_data *kernel_w){

    
    double complex vis = 0 + I*0;
    for(int i = 0; i<pn; ++i){
	double l = *(points + i*2);
	double m = *(points + i*2 + 1);
	double n = sqrt(1.0-l*l-m*m) - 1;
       	double a = 1.0;
	if (apply_aa) a *= grid_correction[(int)(((du*l+theta/2)/theta)*101)] *
			  grid_correction[(int)(((du*m+theta/2)/theta)*101)];
	if (apply_aa_w) a *= grid_correction_w[(int)(((dw*n+theta/2)/theta)*101)];
	//if (apply_aa_w) a *= kernel_w->data[((kernel_w->size+1)*kernel_w->oversampling)/2 +
	//				    (int)(dw*kernel_w->oversampling)];
	double exponent = -2 * M_PI * (u*l + v*m + w*n);
	vis += (cos(exponent) + I * sin(exponent))/a;
    }
    return vis;
}

/*uint64_t predict_grid(double complex *points,
		      double *uvw,
		      double d_u,
		      double d_v,
		      double d_w,
		      struct sep_kernel_data *kernel_uv,
		      struct sep_kernel_data *kernel_w){


}*/


void generate_random_points(double *points,
			    int pn,
			    double theta){

    for(int i=0; i<pn; ++i){
	double l = (float)(randombytes_uniform(65536))/65536 - 0.5;
	double m = (float)(randombytes_uniform(65536))/65536 - 0.5;
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
    
    double *randompoints = malloc(2 * points*sizeof(double));
    generate_random_points(randompoints,points,theta);

    for(int i=0; i<points; ++i){

	printf("l: %f, m: %f",randompoints[i*2],randompoints[i*2+1]);

    }
    printf("\n\n");
    for(int i=0; i<sepkern->size; ++i){
	printf("%f ",sepkern->data[i*sepkern->oversampling]);
    }
    printf("\n");
    
    double complex vis = predict_dft(randompoints, points, 1.0, 1.0, 0.0);
    printf("\nVis: Re: %f, Im: %f \n", creal(vis), cimag(vis));
    double complex vis_sze = predict_sze(randompoints, points, 1.0, 1.0, 0.0, 0.03, 0.03, theta, 1, 1, sepkern, sepkernw);
    printf("Vis Sze: Re: %f, Im: %f \n", creal(vis_sze), cimag(vis_sze));
}
