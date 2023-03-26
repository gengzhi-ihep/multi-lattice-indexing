/* amorphous material diffraction simulator		-James Holton and Ken Frankel		9-12-15

example:

gcc -O -o nonBragg nonBragg.c -lm

./nonBragg -stol water.stol -distance 250 -density 1 -thickness 1

./nonBragg -stol water.stol -lambda 1 -dispersion 0.1 -dispstep 3 -distance 100  -detsize 100 -pixel 0.1 \
  -hdiv 0.28 -hdivstep 0.02 -vdiv 0.28 -vdivstep 0.02 \
  -flux 1e12 -exposure 1

The ".stol" file should contain:
stol F
where stol is sin(theta)/lambda and F is the structure factor of the amorphous 
material.  The structure factor is defined as the ratio of the scattered 
amplitude from the "object" to that of a single electron.  For example, the 
forward-scattered structure factor of water is 2.57 electrons.
wavelength (lambda) should be provided in Angstrom
sample thickness, detector distance, detsize and pixel size in mm
density is in g/cm^3
molecular weight should be in g/mol
divergence in mrad
dispersion in percent
phi and osc are in degrees (for the header)
fluence is in photons/meter^2 (integrated exposure time)

Please note that unlike nearBragg, this program does not work in the near field,
so detector distances should always be much larger than the crystal size

 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#ifndef NAN
#define NAN strtod("NAN",NULL)
#endif

#define TRUE 1
#define FALSE 0
#define Avogadro 6.02214179e23

/* read in text file into double arrays at provided addresses */
size_t read_text_file(char *filename, size_t nargs, ... );

/* cubic spline interpolation function */
void polint(double *xa, double *ya, double x, double *y);


/* rotate a 3-vector in space applied in order phix,phiy,phiz*/
double *rotate(double *v, double *new, double phix, double phiy, double phiz);
/* rotate a 3-vector about a unit vector axis */
double *rotate_axis(double *v, double *new, double *axis, double phi);

/* vector cross product where vector magnitude is 0th element */
double *cross_product(double *x, double *y, double *z);
/* vector inner product where vector magnitude is 0th element */
double dot_product(double *x, double *y);
/* compute difference between two vectors */
double vector_diff(double *vector, double *origin_vector, double *new_vector);
/* measure magnitude of vector and put it in 0th element */
double magnitude(double *vector);
/* scale the magnitude of a vector */
double vector_scale(double *vector, double *new_vector, double scale);
/* force the magnitude of vector to given value */
double vector_rescale(double *vector, double *new_vector, double magnitude);
/* make a unit vector pointing in same direction and report magnitude (both args can be same vector) */
double unitize(double *vector, double *new_unit_vector);


/* polarization factor from vectors */
double polarization_factor(double kahn_factor, double *incident, double *diffracted, double *axis);


/* generate unit vector in random direction */
float uniform3Ddev(float *dx, float *dy, float *dz, long *idum);
/* random deviate with Poisson distribution */
float poidev(float xm, long *idum);
/* random deviate with Gaussian distribution */
float gaussdev(long *idum);
/* random deviate with Lorentzian distribution */
float lorentzdev(long *idum);
/* random deviate with triangle-shaped distribution */
float triangledev(long *idum);
/* random deviate with exponential distribution (>0) */
float expdev(long *idum);
/* random deviate with uniform distribution */
float ran1(long *idum);

/* callback for qsort */
int compare_float(const void *ptr1,const void *ptr2);
float fmedian(unsigned int n, float arr[]);
float fmedian_with_rejection(unsigned int n, float arr[],float sigma,float *mad,int *final_n);
float fmedian_absolute_deviation(unsigned int n, float arr[], float median_value);
float fmean_with_rejection(unsigned int starting_points, float arr[], float sigma_cutoff, float *final_rmsd, int *final_n);

/* file stuff */
char *matfilename = NULL;
char *hklfilename = NULL;
char *dumpfilename = "Fdump.bin\0";
char *stolfilename = NULL;
char *imginfilename = NULL;
char *stoloutfilename = "output.stol\0";
char *histoutfilename = "output.hist\0";
char *sourcefilename = NULL;
char *floatfilename = "floatimage.bin\0";
//char *sinfilename = "sinimage.bin\0";
//char *cosfilename = "cosimage.bin\0";
char *intfilename = "intimage.img\0";
char *pgmfilename = "image.pgm\0";
char *noisefilename = "noiseimage.img\0";
FILE *infile = NULL;
FILE *outfile = NULL;
FILE *stoloutfile = NULL;

typedef enum { SAMPLE, BEAM } pivot;

/* frame handling routines */
typedef struct _SMVinfo
{
	char *filename;
	FILE *handle;
	int swap_bytes;
	int header_size;
	int width;
	int height;
	char *header;
	unsigned short int *mmapdata;
} SMVinfo;

/* SMV image handling routines */
SMVinfo GetFrame(char *filename);
double ValueOf( const char *keyword, SMVinfo smvfile);
unsigned char *read_pgm5_bytes(char *filename,unsigned int *returned_width,unsigned int *returned_height);



int main(int argc, char** argv)
{
    /* progress meter stuff */
    long progress_pixel,progress_pixels;
    int progress_meter=1;
    int babble=1;        
    int printout = 0;
    int printout_ypixel,printout_xpixel=-1;
//    int accumulate = 0;

    /* x-ray beam properties */
    double beam_vector[4]  = {0,1,0,0};
    int coherent = 0;
    int far_source = 0;
    int round_div = 1;
    double lambda,*lambda_of,dispersion=0.0,dispstep=-1,lambda0 = 1.0e-10;
    double *source_X,*source_Y,*source_Z,*source_I,*source_lambda;
    double weight;
    int source,sources;
    double source_path,source_distance = 10.0;
    double hdiv,hdivstep=-1.0,hdivrange= -1.0;
    double vdiv,vdivstep=-1.0,vdivrange= -1.0;
    int divsteps=-1,hdivsteps=-1,vdivsteps=-1,dispsteps=-1;
    int hdiv_tic,vdiv_tic,disp_tic;


    /* Thomson cross section (m^2) */
    double r_e_sqr = 7.94079248018965e-30;
    /* incident x-ray fluence in photons/m^2 */
    double fluence = 125932015286227086360700780544.0;
    double flux=0.0,exposure=1.0,beamsize=1e-4;

    /* things needed to calculate the number of molecules */
    double sample_x   = 1e-4;		/* m */
    double sample_y   = 1e-4;		/* m */
    double sample_z   = 1e-4;		/* m */
    double density    = 1.0e6;		/* g/m^3 */
    double molecular_weight = 18.0;	/* g/mol */
    double volume=0.0,molecules = 0.0;
    /* scale factor = F^2*r_e_sqr*fluence*Avogadro*volume*density/molecular_weight 
                           m^2     ph/m^2  /mol      m^3   g/m^3    g/mol   */

    /* detector stuff */
    double pixel_size = 0.1e-3;
    double pixel_pos[4];
    int xpixel,ypixel,xpixels=0,ypixels=0,pixels;
    double distance = 100.0e-3;
    double detsize_x = 102.4e-3;
    double detsize_y = 102.4e-3;
    double xdet_vector[4]  = {0,0,0,1};
    double ydet_vector[4]  = {0,0,-1,0};
    double zdet_vector[4]  = {0,1,0,0};
    double pix0_vector[4]  = {0,0,0,0};
    double detector_rotx=0.0,detector_roty=0.0,detector_rotz=0.0;
    double twotheta_axis[4] = {0,0,1,0};
    pivot detector_pivot = BEAM;
    double detector_twotheta = 0.0;
    double airpath,omega_pixel,omega_Rsqr_pixel,omega_sum;
    int curved_detector = 0;
    int point_pixel= 0;
    double Xbeam=NAN,Ybeam=NAN;
    double Xdet,Ydet,Rdet;
    double Xdet0,Ydet0;
    double Xclose=NAN,Yclose=NAN,close_distance=NAN;
    double ORGX=NAN,ORGY=NAN;
    double adc_offset = 40.0;


    /* scattering vectors */
    double incident[4];
    double diffracted[4],diffracted0[4];
    double scattering[4];
    double stol,twotheta,theta;
    
    /* polarization stuff */
    double polar_vector[4] = {0,0,0,1};
    double vert_vector[4];
    double polar=1.0,polarization=0.0;
    int nopolar = 0;

    /* sampling */
    int steps;
    int roi_xmin=-1,roi_xmax=-1,roi_ymin=-1,roi_ymax=-1;
    int oversample = -1,subx,suby;
    double subpixel_size;

    /* spindle */
    double phi,phi0=0.0,phistep=-1.0,osc=-1.0;
    int phi_tic,phisteps=-1;
    double spindle_vector[4] = {0,0,0,1};

    /* structure factor representation */
    double F,F_bg,*stol_of,*F_of;
    double F_highangle,F_lowangle;
    int stols,nearest=0;
    double stol_file_mult=1.0e10;
    double denom;

    /* intensity stats */
    double I,I_bg,max_I = 0.0;
    double max_I_x = 0.0,max_I_y = 0.0;
    double intfile_scale = 0.0;
    double pgm_scale = 0.0;
    double sum,sumsqr,avg,rms,rmsd;
    int overloads = 0;        

    /* image file data */
    float *floatimage;
//    float *sinimage;
//    float *cosimage;
    unsigned short int *intimage;
    unsigned char *pgmimage;
    SMVinfo imginfile;
    float *imginfileimage;
    float *diffimage;
    float *stolimage;
    float *Fimage,pixel_F;
    int overflows=0;
    int underflows=0;
    int ignore_values=0;
    unsigned short int ignore_value[70000];
    unsigned short int *invalid_pixel;
    int valid_pixels;

    /* median filter stuff */
    unsigned int bin,*pixels_in,*bin_of;
    float **bin_start;
    float median,mad,deviate,sign;
    float sum_arej,avg_arej,sumd_arej,rms_arej,rmsd_arej;

    /* misc variables */
    int i,j,k,n;
    double X,Y,Z,ratio,r;
    double X0,Y0,Z0,d_r;
    double RTD=180.0*M_1_PI;
    double test;
    double vector[4];
    double newvector[4];
        
    long seed;        
    seed = -time((time_t *)0);
//    printf("random number seed = %u\n",seed);
      
    /* special options */
    int calculate_noise = 1;
    int output_pgm = 1;
    int reject_outliers = 0;


    /* check argument list */
    for(i=1; i<argc; ++i)
    {
        if(argv[i][0] == '-')
        {
            /* option specified */
            if(strstr(argv[i], "-img") && (argc >= (i+1)))
            {
                imginfilename = argv[i+1];
            }
        }
    }

 

    /* read in any provided img file */
    if(imginfilename != NULL)
    {
        /* frame handling routines */
	imginfile = GetFrame(imginfilename);
	if(imginfile.header_size > 0) {
	    xpixels = imginfile.width;
	    ypixels = imginfile.height;
	    pixels = xpixels*ypixels;
	    test = ValueOf("PIXEL_SIZE",imginfile);
	    if(! isnan(test)) pixel_size = test/1000.0;
	    detsize_x = pixel_size*xpixels;
	    detsize_y = pixel_size*ypixels;
	    test = ValueOf("DISTANCE",imginfile);
	    if(! isnan(test)) distance = test/1000.0;
	    test = ValueOf("CLOSE_DISTANCE",imginfile);
	    if(! isnan(test)) close_distance = test/1000.0;
	    test = ValueOf("WAVELENGTH",imginfile);
	    if(! isnan(test)) lambda0 = test/1e10;
	    test = ValueOf("BEAM_CENTER_X",imginfile);
	    if(! isnan(test)) Xbeam = test/1000.0;
	    test = ValueOf("BEAM_CENTER_Y",imginfile);
	    if(! isnan(test)) Ybeam = detsize_y - test/1000.0;
	    test = ValueOf("ORGX",imginfile);
	    if(! isnan(test)) ORGX = test;
	    test = ValueOf("ORGY",imginfile);
	    if(! isnan(test)) ORGY = test;
	    test = ValueOf("PHI",imginfile);
	    if(! isnan(test)) phi0 = test/RTD;
	    test = ValueOf("OSC_RANGE",imginfile);
	    if(! isnan(test)) osc = test/RTD;
	    test = ValueOf("TWOTHETA",imginfile);
	    if(! isnan(test)) twotheta = test/RTD;
	
	    imginfileimage = calloc(pixels+10,sizeof(float));
	    diffimage = calloc(2*pixels+10,sizeof(float));
            stolimage = calloc(pixels+10,sizeof(float));
            Fimage = calloc(pixels+10,sizeof(float));

	    j = imginfile.header_size / sizeof(unsigned short int);
            for(i=0;i<pixels;++i){
	        imginfileimage[i] = (float) imginfile.mmapdata[j];
	         ++j;
	    }
	}
    }

     
    /* check argument list for options */
    for(i=1; i<argc; ++i)
    {
        if(argv[i][0] == '-')
        {
            /* option specified */
            if((strstr(argv[i], "-samplesize") || strstr(argv[i], "-sample_size")) && (argc >= (i+1)))
            {
                sample_x = atof(argv[i+1])/1000;
                sample_y = atof(argv[i+1])/1000;
                sample_z = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-sample_thick") || strstr(argv[i], "-sample_x") || strstr(argv[i], "-thick")) && (argc >= (i+1)))
            {
                sample_x = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-sample_width") || strstr(argv[i], "-sample_y")  || strstr(argv[i], "-width")) && (argc >= (i+1)))
            {
                sample_y = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-sample_heigh") || strstr(argv[i], "-sample_z")  || strstr(argv[i], "-heigh")) && (argc >= (i+1)))
            {
                sample_z = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-density") || strstr(argv[i], "-sample_den")) && (argc >= (i+1)))
            {
                density = atof(argv[i+1])*1e6;
            }
            if((strstr(argv[i], "-molecules") || strstr(argv[i], "-sample_molecules")) && (argc >= (i+1)))
            {
                molecules = atof(argv[i+1]);
            }
            if((0==strcmp(argv[i], "-MW") || strstr(argv[i], "-molecular")) && (argc >= (i+1)))
            {
                molecular_weight = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-Xbeam") && (argc >= (i+1)))
            {
                Xbeam = atof(argv[i+1])/1000.0;
		detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-Ybeam") && (argc >= (i+1)))
            {
                Ybeam = atof(argv[i+1])/1000.0;
		detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-Xclose") && (argc >= (i+1)))
            {
                Xclose = atof(argv[i+1])/1000.0;
		detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-Yclose") && (argc >= (i+1)))
            {
                Yclose = atof(argv[i+1])/1000.0;
		detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-ORGX") && (argc >= (i+1)))
            {
                ORGX = atof(argv[i+1]);
		detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-ORGY") && (argc >= (i+1)))
            {
                ORGY = atof(argv[i+1]);
		detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-pivot") && (argc >= (i+1)))
            {
                if(strstr(argv[i+1], "sample")) detector_pivot = SAMPLE;
		if(strstr(argv[i+1], "beam")) detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-xdet_vector") && (argc >= (i+3)))
            {
                xdet_vector[1] = atof(argv[i+1]);
                xdet_vector[2] = atof(argv[i+2]);
                xdet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-ydet_vector") && (argc >= (i+3)))
            {
                ydet_vector[1] = atof(argv[i+1]);
                ydet_vector[2] = atof(argv[i+2]);
                ydet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-zdet_vector") && (argc >= (i+3)))
            {
                zdet_vector[1] = atof(argv[i+1]);
                zdet_vector[2] = atof(argv[i+2]);
                zdet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-beam_vector") && (argc >= (i+3)))
            {
                beam_vector[1] = atof(argv[i+1]);
                beam_vector[2] = atof(argv[i+2]);
                beam_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-polar_vector") && (argc >= (i+3)))
            {
                polar_vector[1] = atof(argv[i+1]);
                polar_vector[2] = atof(argv[i+2]);
                polar_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-spindle_axis") && (argc >= (i+3)))
            {
                spindle_vector[1] = atof(argv[i+1]);
                spindle_vector[2] = atof(argv[i+2]);
                spindle_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-twotheta_axis") && (argc >= (i+3)))
            {
                twotheta_axis[1] = atof(argv[i+1]);
                twotheta_axis[2] = atof(argv[i+2]);
                twotheta_axis[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-pix0_vector") && (argc >= (i+3)))
            {
                pix0_vector[0] = 1.0;
                pix0_vector[1] = atof(argv[i+1]);
                pix0_vector[2] = atof(argv[i+2]);
                pix0_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-distance") && (argc >= (i+1)))
            {
                distance = atof(argv[i+1])/1000.0;
		detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-close_distance") && (argc >= (i+1)))
            {
                close_distance = atof(argv[i+1])/1000.0;
		detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-source_distance") && (argc >= (i+1)))
            {
		source_distance = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-twotheta") && (argc >= (i+1)))
            {
                detector_twotheta = atof(argv[i+1])/RTD;
		detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-detector_rotx") && (argc >= (i+1)))
            {
                detector_rotx = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-detector_roty") && (argc >= (i+1)))
            {
                detector_roty = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-detector_rotz") && (argc >= (i+1)))
            {
                detector_rotz = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-detsize") && (strlen(argv[i]) == 8) && (argc >= (i+1)))
            {
                detsize_x = atof(argv[i+1])/1000.0;
                detsize_y = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detsize_x") && (argc >= (i+1)))
            {
                detsize_x = atof(argv[i+1])/1000.0;
            }
             if(strstr(argv[i], "-detsize_y") && (argc >= (i+1)))
            {
                detsize_y = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detpixels") && (strlen(argv[i]) == 10) && (argc >= (i+1)))
            {
                xpixels = ypixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detpixels_x") && (argc >= (i+1)))
            {
                xpixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detpixels_y") && (argc >= (i+1)))
            {
                ypixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-curved_det") && (argc >= (i+1)))
            {
                curved_detector = 1;
            }
            if(strstr(argv[i], "-pixel") && (argc >= (i+1)))
            {
                pixel_size = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-point_pixel") )
            {
                point_pixel = 1;
            }
            if(strstr(argv[i], "-polar") && (strlen(argv[i]) == 6) && (argc >= (i+1)))
            {
                polarization = atof(argv[i+1]);
		nopolar = 0;
            }
            if(strstr(argv[i], "-nopolar") )
            {
                nopolar = 1;
            }
            if(strstr(argv[i], "-oversample") && (argc >= (i+1)))
            {
                oversample = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-roi") && (argc >= (i+4)))
            {
                roi_xmin = atoi(argv[i+1]);
                roi_xmax = atoi(argv[i+2]);
                roi_ymin = atoi(argv[i+3]);
                roi_ymax = atoi(argv[i+4]);
            }
            if((strstr(argv[i], "-lambda") || strstr(argv[i], "-wave")) && (argc >= (i+1)))
            {
                lambda0 = atof(argv[i+1])/1.0e10;
            }
            if(strstr(argv[i], "-energy") && (argc >= (i+1)))
            {
                lambda0 = (12398.42/atof(argv[i+1]))/1.0e10;
            }
            if(strstr(argv[i], "-fluence") && (argc >= (i+1)))
            {
                fluence = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-flux") && (argc >= (i+1)))
            {
                flux = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-exposure") && (argc >= (i+1)))
            {
                exposure = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-beamsize") && (argc >= (i+1)))
            {
                beamsize = atof(argv[i+1])/1000;
            }
            if(strstr(argv[i], "-dispersion") && (argc >= (i+1)))
            {
                dispersion = atof(argv[i+1])/100.0;
            }
            if(strstr(argv[i], "-dispsteps") && (argc >= (i+1)))
            {
                dispsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-divergence") && (argc >= (i+1)))
            {
                hdivrange = vdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivrange") && (argc >= (i+1)))
            {
                hdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-vdivrange") && (argc >= (i+1)))
            {
                vdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivstep") && (strlen(argv[i]) == 9) && (argc >= (i+1)))
            {
                hdivstep = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivsteps") && (argc >= (i+1)))
            {
                hdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-vdivstep") && (strlen(argv[i]) == 9) && (argc >= (i+1)))
            {
                vdivstep = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-vdivsteps") && (argc >= (i+1)))
            {
                vdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-divsteps") && (argc >= (i+1)))
            {
                hdivsteps = vdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-round_div") )
            {
		/* cut to circle */
                round_div = 1;
            }
            if(strstr(argv[i], "-square_div") )
            {
		/* just raster */
                round_div = 0;
            }
            if(strstr(argv[i], "-adc") && (argc >= (i+1)))
            {
                adc_offset = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-phi") && strlen(argv[i])==4 && (argc >= (i+1)))
            {
                phi0 = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-osc") && (argc >= (i+1)))
            {
                osc = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-phistep") && strlen(argv[i])==8 && (argc >= (i+1)))
            {
                phistep = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-phisteps") && (argc >= (i+1)))
            {
                phisteps = atoi(argv[i+1]);
            }
//            if(strstr(argv[i], "-dmin") && (argc >= (i+1)))
//            {
//                dmin = atof(argv[i+1])*1e-10;
//            }
//            if(strstr(argv[i], "-mat") && (argc >= (i+1)))
//            {
//                matfilename = argv[i+1];
//            }
//            if(strstr(argv[i], "-hkl") && (argc >= (i+1)))
//            {
//                hklfilename = argv[i+1];
//            }
            if(strstr(argv[i], "-img") && (argc >= (i+1)))
            {
                imginfilename = argv[i+1];
            }
            if(strstr(argv[i], "-ignore") && (argc >= (i+1)))
            {
                ++ignore_values;
                ignore_value[ignore_values] = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-stolout") && strlen(argv[i])>7 && (argc >= (i+1)))
            {
                stoloutfilename = argv[i+1];
            }
            if(strstr(argv[i], "-stol") && strlen(argv[i])==5 && (argc >= (i+1)))
            {
                stolfilename = argv[i+1];
		stol_file_mult = 1e10;
            }
            if(strstr(argv[i], "-4stol") && strlen(argv[i])==6 && (argc >= (i+1)))
            {
                stolfilename = argv[i+1];
		stol_file_mult = 1e10/4;
            }
            if(strstr(argv[i], "-Q") && strlen(argv[i])==2 && (argc >= (i+1)))
            {
                stolfilename = argv[i+1];
		stol_file_mult = 1e10/M_PI/4;
            }
            if(strstr(argv[i], "-sourcefile") && (argc >= (i+1)))
            {
                sourcefilename = argv[i+1];
            }
            if((strstr(argv[i], "-floatfile") || strstr(argv[i], "-floatimage")) && (argc >= (i+1)))
            {
                floatfilename = argv[i+1];
            }
            if((strstr(argv[i], "-intfile") || strstr(argv[i], "-intimage")) && (argc >= (i+1)))
            {
                intfilename = argv[i+1];
            }
            if((strstr(argv[i], "-pgmfile") || strstr(argv[i], "-pgmimage")) && (argc >= (i+1)))
            {
                pgmfilename = argv[i+1];
            }
            if((strstr(argv[i], "-noisefile") || strstr(argv[i], "-noiseimage")) && (argc >= (i+1)))
            {
                noisefilename = argv[i+1];
		calculate_noise = 1;
            }
            if(strstr(argv[i], "-nonoise") )
            {
                /* turn off noise */
                calculate_noise = 0;
            }
            if(strstr(argv[i], "-nopgm") )
            {
                /* turn off noise */
                output_pgm = 0;
            }
            if(strstr(argv[i], "-noreject") )
            {
                /* turn off outlier rejection */
                reject_outliers = 0;
            }
            if(strstr(argv[i], "-scale") && (argc >= (i+1)))
            {
		/* specify the scale for the intfile */
                intfile_scale = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-pgmscale") && (argc >= (i+1)))
            {
		/* specify the scale for the intfile */
                pgm_scale = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-coherent") )
            {
		/* turn off incoherent addition */
                coherent = 1;
            }
            if(strstr(argv[i], "-printout") )
            {
		/* turn on console printing */
                printout = 1;
            }
            if(strstr(argv[i], "-noprogress") )
            {
		/* turn off progress meter */
                progress_meter = 0;
            }
            if(strstr(argv[i], "-progress") )
            {
		/* turn off progress meter */
                progress_meter = 1;
            }
            if(strstr(argv[i], "-printout_pixel") && (argc >= (i+2)))
            {
                printout_xpixel = atoi(argv[i+1]);
                printout_ypixel = atoi(argv[i+2]);
            }
            if(strstr(argv[i], "-seed") && (argc >= (i+1)))
            {
                seed = -atoi(argv[i+1]);
            }
        }
    }

    printf("nonBragg amorphous material diffraction simulator - James Holton and Ken Frankel 3-20-15\n");

    if(stolfilename == NULL){
	printf("usage: nonBragg -stol water.stol\n");
	printf("options:\n");\
	printf("\t-stol filename.stol\ttext file containing sin(theta)/lambda and F for one molecule\n");
        printf("\t-thickness\tthickness of the sample in mm\n");
        printf("\t-samplesize\tlinear dimension of the (cube shaped) sample in mm\n");
        printf("\t-density\tdensity of the sample in g/cm^3\n");
        printf("\t-MW\tmolecular weight of the sample material in g/mol\n");
        printf("\t-hdivrange\thorizontal angular spread of source points in mrad\n");
        printf("\t-vdivrange\tvertical angular spread of source points in mrad\n");
        printf("\t-hdivstep\tnumber of source points in the horizontal\n");
        printf("\t-vdivstep\tnumber of source points in the vertical\n");
        printf("\t-distance\tdistance from origin to detector center in mm\n");
        printf("\t-detsize\tdetector size in mm\n");
        printf("\t-pixel\tdetector pixel size in mm\n");
        printf("\t-oversample\tnumber of sub-pixels per pixel\n");
        printf("\t-lambda\tincident x-ray wavelength in Angstrom\n");
        printf("\t-dispersion\tspectral dispersion: delta-lambda/lambda in percent\n");
        printf("\t-dispsteps\tnumber of wavelengths in above range\n");
        printf("\t-fluence\tintegrated beam intensity in photons/m^2\n");
        printf("\t-flux\t beam intensity in photons/s\n");
        printf("\t-exposure\t exposure time in s\n");
        printf("\t-beamsize\t linear dimension of the (square) beam profile in mm\n");
	printf("\t-sourcefile filename.txt\ttext file containing source positions in mm\n");
        printf("\t-floatfile\tname of binary output file (4-byte floats)\n");
        printf("\t-intfile\tname of smv-formatted output file (arbitrary scale)\n");
        printf("\t-pgmfile\tname of pgm-formatted output file (arbitrary scale)\n");
        printf("\t-noisefile\tname of smv-formatted output file (with Poisson noise)\n");
        printf("\t-Xbeam\timage X coordinate of direct-beam spot (mm)\n");
        printf("\t-Ybeam\timage Y coordinate of direct-beam spot (mm)\n");
        printf("\t-printout\tprint pixel values out to the screen\n");
        printf("\t-noprogress\tturn off the progress meter\n");
	exit(9);
    }


    /* allocate detector memory */
    if(xpixels) {
	detsize_x = pixel_size*xpixels;
    }
    if(ypixels) {
	detsize_y = pixel_size*ypixels;
    }
    xpixels = ceil(detsize_x/pixel_size-0.5);
    ypixels = ceil(detsize_y/pixel_size-0.5);
    pixels = xpixels*ypixels;
    floatimage = calloc(pixels+10,sizeof(float));
    //sinimage = calloc(pixels+10,2*sizeof(float));
    //cosimage = calloc(pixels+10,2*sizeof(float));
    invalid_pixel = calloc(pixels+10,sizeof(unsigned short int));
    intimage   = calloc(pixels+10,sizeof(unsigned short int));
    pgmimage   = calloc(pixels+10,sizeof(unsigned char));

    /* defaults? */
    if(! isnan(ORGX)) Xclose = ORGX*pixel_size;
    if(! isnan(ORGY)) Yclose = ORGY*pixel_size;
    if(isnan(Xclose)) Xclose = (detsize_x - pixel_size)/2.0;
    if(isnan(Yclose)) Yclose = (detsize_y + pixel_size)/2.0;
    if(isnan(Xbeam)) Xbeam = Xclose;
    if(isnan(Ybeam)) Ybeam = Yclose;
    if(roi_xmin < 0) roi_xmin = 0;
    if(roi_xmax < 0) roi_xmax = xpixels;
    if(roi_ymin < 0) roi_ymin = 0;
    if(roi_ymax < 0) roi_ymax = ypixels;
    progress_pixels = (roi_xmax-roi_xmin+1)*(roi_ymax-roi_ymin+1);

    
    if(flux != 0.0 && exposure > 0.0 && beamsize >= 0){
    	fluence = flux*exposure/beamsize/beamsize;
    }
    if(beamsize >= 0){
    	if(beamsize < sample_y){
    	    printf("WARNING: clipping sample (%lg m high) with beam (%lg m)\n",sample_y,beamsize);
    	    sample_y = beamsize;
	}
    	if(beamsize < sample_z){
    	    printf("WARNING: clipping sample (%lg m wide) with beam (%lg m)\n",sample_z,beamsize);
    	    sample_z = beamsize;
	}
    }
    
    /* straighten up sample properties */
    volume = sample_x*sample_y*sample_z;
    if(molecules!=0)
    {
	density = molecules/volume/Avogadro*molecular_weight;
    }
    molecules = volume*density*Avogadro/molecular_weight;

    /* straighten up vectors */
    unitize(beam_vector,beam_vector);
    unitize(xdet_vector,xdet_vector);
    unitize(ydet_vector,ydet_vector);
    cross_product(xdet_vector,ydet_vector,zdet_vector);
    unitize(zdet_vector,zdet_vector);
    unitize(polar_vector,polar_vector);
    unitize(spindle_vector,spindle_vector);
    cross_product(beam_vector,polar_vector,vert_vector);
    unitize(vert_vector,vert_vector);

    
    /* default sampling logic */
    if(phisteps < 0){
        /* auto-select number of phi steps */
	if(osc < 0.0) {
            /* auto-select osc range */
            if(phistep <= 0.0) {
	        /* user doesn't care about anything */
	        phisteps = 1;
		osc = 0.0;
		phistep = 0.0;
	    } else {
		/* user doesn't care about osc or steps, but specified step */
		osc = phistep;
		phisteps = 2;
	    }
	} else {
	    /* user-speficied oscillation */
	    if(phistep <= 0.0) {
	        /* osc specified, but nothing else */
                phisteps = 2;
                phistep = osc/2.0;
            } else {
                /* osc and phi step specified */
                phisteps = ceil(osc/phistep);
	    }
	}
    } else {
	/* user-specified number of phi steps */
	if(phisteps == 0) phisteps = 1;
	if(osc < 0.0) {
            /* auto-select osc range */
            if(phistep <= 0.0) {
	        /* user cares only about number of steps */
		osc = 1.0/RTD;
		phistep = osc/phisteps;
	    } else {
		/* user doesn't care about osc, but specified step */
		osc = phistep;
		phisteps = 2;
	    }
	} else {
	    /* user-speficied oscillation */
	    if(phistep < 0.0) {
	        /* osc and steps specified */
                phistep = osc/phisteps;
            } else {
                /* everything specified */
	    }
	}
    }

    if(hdivsteps <= 0){
        /* auto-select number of steps */
	if(hdivrange < 0.0) {
            /* auto-select range */
            if(hdivstep <= 0.0) {
	        /* user doesn't care about anything */
	        hdivsteps = 1;
		hdivrange = 0.0;
		hdivstep = 0.0;
	    } else {
		/* user specified stepsize and nothing else */
		hdivrange = hdivstep;
		hdivsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(hdivstep <= 0.0) {
	        /* range specified, but nothing else */
                hdivstep = hdivrange;
                hdivsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                hdivsteps = ceil(hdivrange/hdivstep);
	    }
	}
    } else {
	/* user-specified number of steps */
	if(hdivrange < 0.0) {
            /* auto-select range */
            if(hdivstep <= 0.0) {
	        /* user cares only about number of steps */
		hdivrange = 1.0;
		hdivstep = hdivrange/hdivsteps;
	    } else {
		/* user doesn't care about range */
		hdivrange = hdivstep;
		hdivsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(hdivstep <= 0.0) {
	        /* range and steps specified */
		if(hdivsteps <=1 ) hdivsteps = 2;
                hdivstep = hdivrange/(hdivsteps-1);
            } else {
                /* everything specified */
	    }
	}
    }

    if(vdivsteps <= 0){
        /* auto-select number of steps */
	if(vdivrange < 0.0) {
            /* auto-select range */
            if(vdivstep <= 0.0) {
	        /* user doesn't care about anything */
	        vdivsteps = 1;
		vdivrange = 0.0;
		vdivstep = 0.0;
	    } else {
		/* user specified stepsize and nothing else */
		vdivrange = vdivstep;
		vdivsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(vdivstep <= 0.0) {
	        /* range specified, but nothing else */
                vdivstep = vdivrange;
                vdivsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                vdivsteps = ceil(vdivrange/vdivstep);
	    }
	}
    } else {
	/* user-specified number of steps */
	if(vdivrange < 0.0) {
            /* auto-select range */
            if(vdivstep <= 0.0) {
	        /* user cares only about number of steps */
		vdivrange = 1.0;
		vdivstep = vdivrange/vdivsteps;
	    } else {
		/* user doesn't care about range */
		vdivrange = vdivstep;
		vdivsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(vdivstep <= 0.0) {
	        /* range and steps specified */
		if(vdivsteps <=1 ) vdivsteps = 2;
                vdivstep = vdivrange/(vdivsteps-1);
            } else {
                /* everything specified */
	    }
	}
    }
    

    if(dispsteps <= 0){
        /* auto-select number of steps */
	if(dispersion < 0.0) {
            /* auto-select range */
            if(dispstep <= 0.0) {
	        /* user doesn't care about anything */
	        dispsteps = 1;
		dispersion = 0.0;
		dispstep = 0.0;
	    } else {
		/* user specified stepsize and nothing else */
		dispersion = dispstep;
		dispsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(dispstep <= 0.0) {
	        /* range specified, but nothing else */
                dispstep = dispersion;
                dispsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                dispsteps = ceil(dispersion/dispstep);
	    }
	}
    } else {
	/* user-specified number of steps */
	if(dispersion < 0.0) {
            /* auto-select range */
            if(dispstep <= 0.0) {
	        /* user cares only about number of steps */
		dispersion = 1.0;
		dispstep = dispersion/dispsteps;
	    } else {
		/* user doesn't care about range */
		dispersion = dispstep;
		dispsteps = 2;
	    }
	} else {
	    /* user-speficied range */
	    if(dispstep <= 0.0) {
	        /* range and steps specified */
		if(dispsteps <=1 ) dispsteps = 2;
                dispstep = dispersion/(dispsteps-1);
            } else {
                /* everything specified */
	    }
	}
    }
        
    
    /* sanity checks */
    if(hdivrange <= 0.0 || hdivstep <= 0.0 || hdivsteps <= 0) {
        hdivsteps = 1;
        hdivrange = 0.0;
        hdivstep = 0.0;
    }
    if(vdivrange <= 0.0 || vdivstep <= 0.0 || vdivsteps <= 0) {
        vdivsteps = 1;
        vdivrange = 0.0;
        vdivstep = 0.0;
    }
    if(dispersion <= 0.0 || dispstep <= 0.0 || dispsteps <= 0) {
        dispsteps = 1;
        dispersion = 0.0;
        dispstep = 0.0;
    }

   
    
    /* initialize detector origin from a beam center and distance */
    /* there are two conventions here: mosflm and XDS */

    /* first off, what is the relationship between the two "beam centers"? */
    rotate(zdet_vector,vector,detector_rotx,detector_roty,detector_rotz);
    ratio = dot_product(beam_vector,vector);
    if(ratio == 0.0) { ratio = DBL_MIN; }
    if(isnan(close_distance)) close_distance = ratio*distance;
    distance = close_distance/ratio;

    if(detector_pivot == SAMPLE){
        printf("pivoting detector around sample\n");
	/* initialize detector origin before rotating detector */
        pix0_vector[1] = -Xclose*xdet_vector[1]-Yclose*ydet_vector[1]+close_distance*zdet_vector[1];
        pix0_vector[2] = -Xclose*xdet_vector[2]-Yclose*ydet_vector[2]+close_distance*zdet_vector[2];
        pix0_vector[3] = -Xclose*xdet_vector[3]-Yclose*ydet_vector[3]+close_distance*zdet_vector[3];

        /* now swing the detector origin around */
        rotate(pix0_vector,pix0_vector,detector_rotx,detector_roty,detector_rotz);
        rotate_axis(pix0_vector,pix0_vector,twotheta_axis,detector_twotheta);
    }
    /* now orient the detector plane */
    rotate(xdet_vector,xdet_vector,detector_rotx,detector_roty,detector_rotz);
    rotate(ydet_vector,ydet_vector,detector_rotx,detector_roty,detector_rotz);
    rotate(zdet_vector,zdet_vector,detector_rotx,detector_roty,detector_rotz);

    /* also apply orientation part of twotheta swing */
    rotate_axis(xdet_vector,xdet_vector,twotheta_axis,detector_twotheta);
    rotate_axis(ydet_vector,ydet_vector,twotheta_axis,detector_twotheta);
    rotate_axis(zdet_vector,zdet_vector,twotheta_axis,detector_twotheta);
    
    /* make sure beam center is preserved */
    if(detector_pivot == BEAM){
        printf("pivoting detector around direct beam spot\n");
        pix0_vector[1] = -Xbeam*xdet_vector[1]-Ybeam*ydet_vector[1]+distance*beam_vector[1];
        pix0_vector[2] = -Xbeam*xdet_vector[2]-Ybeam*ydet_vector[2]+distance*beam_vector[2];
        pix0_vector[3] = -Xbeam*xdet_vector[3]-Ybeam*ydet_vector[3]+distance*beam_vector[3];
    }

    /* what is the point of closest approach between sample and detector? */
    Xclose         = -dot_product(pix0_vector,xdet_vector);
    Yclose         = -dot_product(pix0_vector,ydet_vector);
    close_distance =  dot_product(pix0_vector,zdet_vector);

    /* where is the direct beam now? */
    /* difference between beam impact vector and detector origin */
    newvector[1] = close_distance/ratio*beam_vector[1]-pix0_vector[1];
    newvector[2] = close_distance/ratio*beam_vector[2]-pix0_vector[2];
    newvector[3] = close_distance/ratio*beam_vector[3]-pix0_vector[3];
    /* extract components along detector vectors */
    Xbeam = dot_product(xdet_vector,newvector);
    Ybeam = dot_product(ydet_vector,newvector);
    distance = close_distance/ratio;    

        

    /* now read in amorphous material structure factors */
    printf("reading %s\n",stolfilename);
    stols = read_text_file(stolfilename,2,&stol_of,&F_of);
    if(stols == 0){
    	perror("no data in input file");
	exit(9);
    }
    /* add two values at either end for interpolation */
    stols += 4;
    F_highangle = NAN;
    for(i=stols-3;i>1;--i){
	stol_of[i] = stol_of[i-2] * stol_file_mult;
	F_of[i]    = F_of[i-2];
	if(! isnan(F_of[i])) {
	    F_lowangle = F_of[i];
	    if(isnan(F_highangle)) {
		F_highangle = F_of[i];
	    }
	}
	else
	{
	    /* missing values are zero */
	    F_of[i] = 0.0;
	}
    }
    stol_of[0] = -1e99;
    stol_of[1] = -1e98;
    F_of[0] = F_of[1] = F_lowangle;
    stol_of[stols-2] = 1e98;
    stol_of[stols-1] = 1e99;
    F_of[stols-1] = F_of[stols-2] = F_highangle;

    /* allocate memory for counting how many of these get used */
    bin_start = calloc(stols,sizeof(float **));
    pixels_in = calloc(stols,sizeof(unsigned int));
    bin_of    = calloc(pixels+10,sizeof(unsigned int));
   
    /* import sources from user file */
    sources = 0;
    if(sourcefilename != NULL) {
        sources = read_text_file(sourcefilename,5,&source_X,&source_Y,&source_Z,&source_I,&source_lambda);
	if(sources == 0) {
	    perror("reading source definition file");
	    exit(9);
	}
        /* apply defaults to missing values */
        for(source=0;source<sources;++source){
            if(isnan(source_X[source])) {
                source_X[source] = -source_distance*beam_vector[1];
            }
            if(isnan(source_Y[source])) {
                source_Y[source] = -source_distance*beam_vector[2];
            }
            if(isnan(source_Z[source])) {
                source_Z[source] = -source_distance*beam_vector[3];
            }
            if(isnan(source_I[source])) {
                source_I[source] = 1.0;
            }
            if(isnan(source_lambda[source])) {
                source_lambda[source] = lambda0;
            }
	}
    }
   
   
    if(sources == 0)
    {
    	/* generate generic list of sources */
    
        /* count divsteps sweep over solid angle of beam divergence */
        divsteps = 0;
        for(hdiv_tic=0;hdiv_tic<hdivsteps;++hdiv_tic){
            for(vdiv_tic=0;vdiv_tic<vdivsteps;++vdiv_tic){
                hdiv = hdivstep * hdiv_tic - hdivrange/2.0 ;
                vdiv = vdivstep * vdiv_tic - vdivrange/2.0 ;
                /* force an elliptical divergence */
                test = (hdiv*hdiv-hdivstep*hdivstep/4.0*(1-hdivsteps%2))/hdivrange/hdivrange ;
                test += (vdiv*vdiv-vdivstep*vdivstep/4.0*(1-vdivsteps%2))/vdivrange/vdivrange ;
                if( round_div && test*4.0 > 1.1) continue;

                ++divsteps;
                printf("divergence deviation: %g %g\n",hdiv,vdiv);
            }
        }

        /* print out wavelength steps with sweep over spectral dispersion */
        for(disp_tic=0;disp_tic<dispsteps;++disp_tic){
            lambda = lambda0 * ( 1.0 + dispstep * disp_tic - dispersion/2.0 ) ;
            printf("lambda%d = %.15g\n",disp_tic,lambda);
        }
	
	/* allocate enough space */
	sources = divsteps*dispsteps;
	source_X = calloc(sources+10,sizeof(double));
	source_Y = calloc(sources+10,sizeof(double));
	source_Z = calloc(sources+10,sizeof(double));
	source_I = calloc(sources+10,sizeof(double));
	source_lambda = calloc(sources+10,sizeof(double));
	
	/* now actually create the source entries */
	sources = 0;
        for(hdiv_tic=0;hdiv_tic<hdivsteps;++hdiv_tic){
            for(vdiv_tic=0;vdiv_tic<vdivsteps;++vdiv_tic){
                hdiv = hdivstep * hdiv_tic - hdivrange/2.0 ;
                vdiv = vdivstep * vdiv_tic - vdivrange/2.0 ;
                /* force an elliptical divergence */
                test = (hdiv*hdiv-hdivstep*hdivstep/4.0*(1-hdivsteps%2))/hdivrange/hdivrange ;
                test += (vdiv*vdiv-vdivstep*vdivstep/4.0*(1-vdivsteps%2))/vdivrange/vdivrange ;
                if( round_div && test*4.0 > 1.1) continue;

	        /* construct unit vector along "beam" */
	        vector[1] = -source_distance*beam_vector[1];
	        vector[2] = -source_distance*beam_vector[2];
	        vector[3] = -source_distance*beam_vector[3];
	        /* divergence is in angle space */
		/* define "horizontal" as the E-vector of the incident beam */
                rotate_axis(vector,newvector,polar_vector,vdiv);
        	rotate_axis(newvector,vector,vert_vector,hdiv);

		/* one source at each position for each wavelength */
	        for(disp_tic=0;disp_tic<dispsteps;++disp_tic){
	            lambda = lambda0 * ( 1.0 + dispstep * disp_tic - dispersion/2.0 ) ;

		    source_X[sources] = vector[1];
		    source_Y[sources] = vector[2];
		    source_Z[sources] = vector[3];
		    source_I[sources] = 1.0;
		    source_lambda[sources] = lambda;
		    ++sources;
		}
    	    }
    	}
    }
    printf("  created a total of %d sources:\n",sources);
    for(source=0;source<sources;++source){

    	/* retrieve stuff from cache */
	X = source_X[source];
	Y = source_Y[source];
	Z = source_Z[source];
	I = source_I[source];
	lambda = source_lambda[source];

    	printf("%g %g %g   %g %g\n",X,Y,Z,I,lambda);
    }


    /* final decisions about sampling */
    if(oversample <= 0) oversample = 1;
    steps = sources*oversample*oversample;
    subpixel_size = pixel_size/oversample;
 

    printf("  %d initialized F points (will cubic-spline interpolate between them)\n",stols);
    printf("  wave=%g meters +/- %g%% in %d steps\n",lambda0,dispersion*100,dispsteps);
    if(nopolar) { printf("  polarization effect disabled\n"); }
           else { printf("  Kahn polarization factor: %f\n",polarization); }
    if(curved_detector) printf("  curved detector: all pixels same distance from origin\n");
    if(point_pixel) printf("  pixel obliquity effect disabled\n");
    printf("  incident fluence: %lg photons/m^2\n",fluence);
    printf("  sample is %lg m thick x %lg m high x %lg m wide, %lg g/cm^3 and %lg g/mol (%lg molecules)\n",
            sample_x,sample_y,sample_z,density/1e6,molecular_weight,molecules);
    printf("  distance=%g detsize=%gx%g  pixel=%g meters (%dx%d pixels)\n",distance,detsize_x,detsize_y,pixel_size,xpixels,ypixels);
    printf("  Xbeam=%g Ybeam=%g\n",Xbeam,Ybeam);
    printf("  detector origin: %g %g %g\n",pix0_vector[1],pix0_vector[2],pix0_vector[3]);
    printf("  DIRECTION_OF_DETECTOR_X-AXIS= %g %g %g\n",xdet_vector[1],xdet_vector[2],xdet_vector[3]);
    printf("  DIRECTION_OF_DETECTOR_Y-AXIS= %g %g %g\n",ydet_vector[1],ydet_vector[2],ydet_vector[3]);
    printf("  INCIDENT_BEAM_DIRECTION= %g %g %g\n",beam_vector[1],beam_vector[2],beam_vector[3]);
    printf("  POLARIZATION_PLANE_NORMAL= %g %g %g\n",polar_vector[1],polar_vector[2],polar_vector[3]);
    printf("  roi: %d < x < %d && %d < y < %d\n",roi_xmin,roi_xmax,roi_ymin,roi_ymax);
    printf("  hdivrange=%g hdivstep=%g  radians\n",hdivrange,hdivstep);
    printf("  vdivrange=%g vdivstep=%g  radians\n",vdivrange,vdivstep);
    printf("  %d divergence steps\n",divsteps);
    printf("  %d sources\n",sources);
    printf("  %d pixel oversample steps\n",oversample);
//    printf("  coherent source: %d\n",coherent);
    if(calculate_noise){
        printf("\n  noise image paramters:\n");
        printf("  seed: %ld\n",seed);
    } 


    /* sweep over detector */   
    sum = sumsqr = 0.0;
    j = 0;
    progress_pixel = 0;
    valid_pixels = 0;
    omega_sum = 0.0;
    for(ypixel=0;ypixel<ypixels;++ypixel){
	for(xpixel=0;xpixel<xpixels;++xpixel){

    	    /* allow for just one part of detector to be rendered */
	    if(xpixel < roi_xmin || xpixel > roi_xmax || ypixel < roi_ymin || ypixel > roi_ymax) {
		++invalid_pixel[j];
		++j; continue;
	    }

	    /* reset photon count for this pixel */
	    I = 0;

	    /* loop over sub-pixels */
	    for(suby=0;suby<oversample;++suby){
		for(subx=0;subx<oversample;++subx){

		    /* absolute mm position on detector (relative to its origin) */
		    Xdet = subpixel_size*(xpixel*oversample + subx ) + subpixel_size/2.0;
		    Ydet = subpixel_size*(ypixel*oversample + suby ) + subpixel_size/2.0;
//		    Xdet = pixel_size*xpixel;
//		    Ydet = pixel_size*ypixel;

		    /* construct detector pixel position in 3D space */
//		    pixel_X = distance;
//		    pixel_Y = Ydet-Ybeam;
//		    pixel_Z = Xdet-Xbeam;
        	    pixel_pos[1] = Xdet*xdet_vector[1]+Ydet*ydet_vector[1]+pix0_vector[1];
        	    pixel_pos[2] = Xdet*xdet_vector[2]+Ydet*ydet_vector[2]+pix0_vector[2];
        	    pixel_pos[3] = Xdet*xdet_vector[3]+Ydet*ydet_vector[3]+pix0_vector[3];
                    pixel_pos[0] = 0.0;
		    if(curved_detector) {
			/* construct detector pixel that is always "distance" from the sample */
			vector[1] = distance*beam_vector[1]; vector[2]=distance*beam_vector[2] ; vector[3]=distance*beam_vector[3];
			/* treat detector pixel coordinates as radians */
        		rotate_axis(vector,newvector,ydet_vector,pixel_pos[2]/distance);
        		rotate_axis(newvector,pixel_pos,xdet_vector,pixel_pos[3]/distance);
// 	    		rotate(vector,pixel_pos,0,pixel_pos[3]/distance,pixel_pos[2]/distance);
		    }
		    /* construct the diffracted-beam unit vector to this pixel */
		    airpath = unitize(pixel_pos,diffracted);

		    /* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
		    omega_pixel = pixel_size*pixel_size/airpath/airpath*close_distance/airpath;
		    /* option to turn off obliquity effect, inverse-square-law only */
                    if(point_pixel) omega_pixel = 1.0/airpath/airpath;
		    omega_sum += omega_pixel;

		    /* loop over sources now */
		    for(source=0;source<sources;++source){

    	    	    	/* retrieve stuff from cache */
			incident[1] = -source_X[source];
			incident[2] = -source_Y[source];
			incident[3] = -source_Z[source];
			lambda = source_lambda[source];

			/* construct the incident beam unit vector while recovering source distance */
			source_path = unitize(incident,incident);

			/* construct the scattering vector for this pixel */
			scattering[1] = (diffracted[1]-incident[1])/lambda;
			scattering[2] = (diffracted[2]-incident[2])/lambda;
			scattering[3] = (diffracted[3]-incident[3])/lambda;

    	    	    	/* sin(theta)/lambda is half the scattering vector length */
			stol = 0.5*magnitude(scattering);

			/* now we need to find the nearest four "stol file" points */
			while(stol > stol_of[nearest] && nearest <= stols){++nearest; };
			while(stol < stol_of[nearest] && nearest >= 2){--nearest; };

			/* cubic spline interpolation */
			polint(stol_of+nearest-1, F_of+nearest-1, stol, &F);
//			if(F<0.0) F=0.0;
			sign=1.0;
			if(F<0.0) sign=-1.0;

    	    	    	/* now we have the structure factor for this pixel */

			/* polarization factor */
			if(! nopolar){
			    /* need to compute polarization factor */
			    polar = polarization_factor(polarization,incident,diffracted,polar_vector);
			}
			else
			{
			    polar = 1.0;
			}

			/* accumulate unscaled pixel intensity from this */
			I += sign*F*F*polar*omega_pixel*source_I[source];
		    }
		    /* end of source loop */
		}
		/* end of sub-pixel y loop */
            }
	    /* end of sub-pixel x loop */


	    /* save photons/pixel (if fluence specified), or F^2/omega if no fluence given */
	    floatimage[j]= I*r_e_sqr*fluence*molecules/steps;
	    
	    if(imginfilename != NULL) {
		/* is the pixel valid on the input image? */
		/* skip over any invalid values */
		for(k=1;k<=ignore_values;++k)
	        {
		    if(imginfileimage[j]==ignore_value[k]){
		        ++invalid_pixel[j];
		    }
	        }
		
		/* transform pixel intensity back to a structure factor */
		deviate=imginfileimage[j]-adc_offset;
		sign = 1.0;
		if(deviate<0.0) sign = -1.0;
		deviate = fabsf(deviate);
	    	pixel_F = sign*sqrt(deviate/polar/omega_pixel/fluence/r_e_sqr/molecules*steps);
		/* maintain F and stol images */
                stolimage[j] = stol/stol_file_mult;
                Fimage[j] = pixel_F;
		bin = 0;
	        if(! invalid_pixel[j])
		{
		    /* figure out which stol bin this pixel belongs to.  invalid pixels are in bin=0 */
		    bin = nearest;
		    if(stol > (stol_of[bin]+stol_of[bin+1])/2.0) ++bin;
		    ++valid_pixels;
		}
		++pixels_in[bin];
		bin_of[j]=bin;
	    }

	    if(floatimage[j] > max_I) {
	        max_I = floatimage[j];
	        max_I_x = Xdet;
	        max_I_y = Ydet;
	    }
	    sum += floatimage[j];
            sumsqr += floatimage[j]*floatimage[j];
            ++n;
	    
	    if( printout )
	    {
		if((xpixel==printout_xpixel && ypixel==printout_ypixel) || printout_xpixel < 0)
		{
		    twotheta = atan2(sqrt(pixel_pos[2]*pixel_pos[2]+pixel_pos[3]*pixel_pos[3]),pixel_pos[1]);
		    test = sin(twotheta/2.0)/(lambda0*1e10);
	    	    printf("%4d %4d : stol = %g or %g\n", xpixel,ypixel,stol,test);
	    	    printf(" F=%g    I = %g\n", F,I);
	    	    printf("I/steps %15.10g\n", I/steps);
	    	    printf("polar   %15.10g\n", polar);
	    	    printf("omega   %15.10g\n", omega_pixel);
	    	    printf("pixel   %15.10g\n", floatimage[j]);
		}
	    }
	    else
	    {
		if(progress_meter && progress_pixels/100 > 0)
		{
	            if(progress_pixel % ( progress_pixels/20 ) == 0 ||
                       ((10*progress_pixel<progress_pixels ||
                         10*progress_pixel>9*progress_pixels) && 
                        (progress_pixel % (progress_pixels/100) == 0)))
		    {
			printf("%lu%% done\n",progress_pixel*100/progress_pixels);
	            }
		}
	    	++progress_pixel;
    	    }
	    ++j;
    	}
    }
    printf("\n");

    printf("solid angle subtended by detector = %g steradian ( %g%% sphere)\n",omega_sum,100*omega_sum/4/M_PI);

    if(imginfilename != NULL && stoloutfilename != NULL)
    {
	outfile = fopen(stoloutfilename,"w");
	if(outfile == NULL) {
	    perror(stoloutfilename);
	    exit(9);
	}
    
	/* set up pointers with enough space after each of them */
	bin_start[0]=calloc(2*pixels+10*stols,sizeof(float));
        ++bin_start[0];
	for(bin=1;bin<stols-1;++bin)
	{
	    /* each array must have 2*n values in it */
	    bin_start[bin]=bin_start[bin-1]+2*pixels_in[bin-1]+2;
	}

	/* populate each bin with appropriate pixel values */
	for(j=0;j<pixels;++j)
	{
	    bin = bin_of[j];
	    *bin_start[bin] = Fimage[j];
	    /* increment the pointer to the next value */
	    ++bin_start[bin];
	    /* we will reset the starting points in the next loop */
	}

        i=0;
	for(bin=2;bin<stols-2;++bin)
	{
	    /* correct pointer drift */
	    bin_start[bin] -= pixels_in[bin];

	    stol = stol_of[bin];
	    /* this function looks at "input_n" elements, starting at 1 */
	    median   = fmedian_with_rejection(pixels_in[bin],bin_start[bin]-1,6.0,&mad,&n);
	    avg_arej = fmean_with_rejection(n,bin_start[bin],6.0,&rmsd_arej,&n);
	    if(n>100)
	    {
//	  	fprintf(outfile,"%g %g %g  %d\n",stol/stol_file_mult,median,mad,n);
		fprintf(outfile,"%g %g %g  %d\n",stol/stol_file_mult,avg_arej,rmsd_arej,n);
		++i;
	    }
	    else
	    {
		printf("WARNING: not enough pixels in bin= %d n=%d stol= %g median= %g avg_arej= %g\n",bin,n,stol/stol_file_mult,median,avg_arej);
	    }
	}
	printf("wrote %s as %d lines of text\n",stoloutfilename,i);
	fclose(outfile);
    }

    /* do some stats? */
    if(n<=0) n=1;
    avg = sum/n;
    if(n<=1) n=2;
    rms = sqrt(sumsqr/(n-1));
    sumsqr = 0.0;
    j = n = 0;
    for(ypixel=0;ypixel<ypixels;++ypixel)
    {
        for(xpixel=0;xpixel<xpixels;++xpixel)
        {
	    ++j;
            if(xpixel < roi_xmin || xpixel > roi_xmax || ypixel < roi_ymin || ypixel > roi_ymax)
            {
		continue;
	    }
	    test = floatimage[j]-avg;
	    sumsqr += test*test;
	    ++n;
        }
    }
    if(n<=1) n=2;
    rmsd = sqrt(sumsqr/(n-1));


    printf("writing %s as %d %lu-byte floats\n",floatfilename,pixels,sizeof(float));
    outfile = fopen(floatfilename,"w");
    if(outfile == NULL)
    {
	perror("ERROR: fopen");
	exit(9);
    }
    fwrite(floatimage,sizeof(float),pixels,outfile);
    fclose(outfile);

    /* output as ints */   
    j = 0;
    printf("max_I = %g  at %g %g\n",max_I,max_I_x,max_I_y);
    printf("mean= %g rms= %g rmsd= %g\n",avg,rms,rmsd);
    if(intfile_scale <= 0.0){
	intfile_scale = 1.0;
	if(max_I > 0.0) intfile_scale = 55000.0/max_I;
    }
    printf("intfile_scale = %g\n",intfile_scale);
    for(ypixel=0;ypixel<ypixels;++ypixel)
    {
        for(xpixel=0;xpixel<xpixels;++xpixel)
        {
            if(xpixel < roi_xmin || xpixel > roi_xmax || ypixel < roi_ymin || ypixel > roi_ymax)
            {
	       ++j; continue;
            }
	    test = floatimage[j] *intfile_scale+adc_offset;
	    if(test > 65535.0) test = 65535.0;
	    if(test < 0.0) test = 0.0;
	    intimage[j] = (unsigned short int) ( floorf(test+0.5) );
//	    printf("%d %d = %d\n",xpixel,ypixel,intimage[j]);
	    ++j;
        }
    }

    printf("writing %s as %lu-byte integers\n",intfilename,sizeof(unsigned short int));
    outfile = fopen(intfilename,"w");
    if(outfile == NULL)
    {
	    perror("ERROR: fopen");
	    exit(9);
    }
    fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=little_endian;\nTYPE=unsigned_short;\n");
    fprintf(outfile,"SIZE1=%d;\nSIZE2=%d;\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",xpixels,ypixels,pixel_size*1000.0,distance*1000.0);
    fprintf(outfile,"WAVELENGTH=%g;\n",lambda0*1e10);
    fprintf(outfile,"BEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",Xbeam*1000.0,(detsize_y-Ybeam)*1000);
    fprintf(outfile,"ORGX=%g;\nORGY=%g;\n",Xclose/pixel_size,Yclose/pixel_size);
    fprintf(outfile,"CLOSE_DISTANCE=%g;\n",close_distance*1000.0);
    fprintf(outfile,"PHI=%g;\nOSC_START=%g;\nOSC_RANGE=%g;\n",phi0*RTD,phi0*RTD,osc*RTD);
    fprintf(outfile,"TWOTHETA=%g;\n",detector_twotheta*RTD);
    fprintf(outfile,"DETECTOR_SN=000;\n");
    fprintf(outfile,"BEAMLINE=fake;\n");
    fprintf(outfile,"}\f");
    while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    fwrite(intimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);

    /* compare to original? */
    valid_pixels = 0;
    if(imginfilename != NULL)
    {
	for(i=0;i<pixels;++i)
	{
	    if(! invalid_pixel[i])
	    {
		++valid_pixels;
		deviate = imginfileimage[i] - floatimage[i];
		diffimage[valid_pixels] = deviate;
	    }
	}
	if(reject_outliers)
	{
	    median   = fmedian_with_rejection(valid_pixels,diffimage-1,6.0,&mad,&n);
	    printf("difference from input image after outlier rejection: median= %g mad= %g ( %d pixels)\n",median,mad,n);
	    avg_arej = fmean_with_rejection(n,diffimage-1,4.0,&rmsd_arej,&n);
	    sumsqr=0.0;
	    for(j=1;j<=n;++j)
	    {
	        sumsqr += diffimage[j]*diffimage[j];
	    }
	    rms_arej=sqrt(sumsqr/n);
	    printf("difference from input image after outlier rejection: mean= %g rms= %g rmsd= %g ( %d pixels)\n",avg_arej,rms_arej,rmsd_arej,n);
	}
    }


    /* output as pgm */   
    j = 0;
    if(pgm_scale <= 0.0){
        pgm_scale = intfile_scale;
	if(rmsd > 0.0) pgm_scale = 250.0/(5.0*rmsd);
    }
    printf("pgm_scale = %g\n",pgm_scale);
    j = 0;
    for(ypixel=0;ypixel<ypixels;++ypixel)
    {
        for(xpixel=0;xpixel<xpixels;++xpixel)
        {
            if(xpixel < roi_xmin || xpixel > roi_xmax || ypixel < roi_ymin || ypixel > roi_ymax)
            {
                ++j; continue;
            }
	    test = floatimage[j] * pgm_scale;
	    if(test > 255.0) test = 255.0;
	    pgmimage[j] = (unsigned char) ( test );
//	    printf("%d %d = %d\n",xpixel,ypixel,pgmimage[j]);
	    ++j;
        }
    }

    printf("writing %s as %lu-byte integers\n",pgmfilename,sizeof(unsigned char));
    outfile = fopen(pgmfilename,"w");
    if(outfile == NULL)
    {
	    perror("ERROR: fopen");
	    exit(9);
    }
    fprintf(outfile, "P5\n%d %d\n", xpixels, ypixels);
    fprintf(outfile, "# pixels scaled by %lg\n", pgm_scale);
    fprintf(outfile, "255\n");
    fwrite(pgmimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);


    if(calculate_noise == 0){
	return 0;
    }

    /* simulate Poisson noise */
    j = 0;
    sum = 0.0;
    overloads = 0;
    for(ypixel=0;ypixel<ypixels;++ypixel)
    {
        for(xpixel=0;xpixel<xpixels;++xpixel)
        {
            if(xpixel < roi_xmin || xpixel > roi_xmax || ypixel < roi_ymin || ypixel > roi_ymax) 
            {
                ++j; continue;
            }
	    test = poidev( floatimage[j], &seed );
	    sum += test;
	    test += adc_offset;
	    if(test > 65535.0)
            {
	        test = 65535.0;
	        ++overloads;
	    }
	    intimage[j] = (unsigned short int) test;
//	    printf("%d %d = %d\n",xpixel,ypixel,intimage[j]);
	    ++j;
        }
    }
    printf("%.0f photons on noise image (%d overloads)\n",sum,overloads);

    printf("writing %s as %lu-byte integers\n",noisefilename,sizeof(unsigned short int));
    outfile = fopen(noisefilename,"w");
    if(outfile == NULL)
    {
	    perror("ERROR: fopen");
	    exit(9);
    }
    fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=little_endian;\nTYPE=unsigned_short;\n");
    fprintf(outfile,"SIZE1=%d;\nSIZE2=%d;\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",xpixels,ypixels,pixel_size*1000.0,distance*1000.0);
    fprintf(outfile,"WAVELENGTH=%g;\n",lambda0*1e10);
    fprintf(outfile,"BEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",Xbeam*1000.0,(detsize_y-Ybeam)*1000);
    fprintf(outfile,"ORGX=%g;\nORGY=%g;\n",Xclose/pixel_size,Yclose/pixel_size);
    fprintf(outfile,"CLOSE_DISTANCE=%g;\n",close_distance*1000.0);
    fprintf(outfile,"PHI=%g;\nOSC_START=%g;\nOSC_RANGE=%g;\n",phi0*RTD,phi0*RTD,osc*RTD);
    fprintf(outfile,"TWOTHETA=%g;\n",detector_twotheta*RTD);
    fprintf(outfile,"DETECTOR_SN=000;\n");
    fprintf(outfile,"BEAMLINE=fake;\n");
    fprintf(outfile,"}\f");
    while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    fwrite(intimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);

    return 0;
}


double *rotate(double *v, double *new, double phix, double phiy, double phiz) {
    
    double rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz;
    double new_x,new_y,new_z,rotated_x,rotated_y,rotated_z;
    
    new_x=v[1];
    new_y=v[2];
    new_z=v[3];
    
    if(phix != 0){
        /* rotate around x axis */
        //rxx= 1;         rxy= 0;         rxz= 0;
        ryx= 0;         ryy= cos(phix); ryz=-sin(phix);
        rzx= 0;         rzy= sin(phix); rzz= cos(phix);
        
        rotated_x = new_x;
        rotated_y = new_y*ryy + new_z*ryz;
        rotated_z = new_y*rzy + new_z*rzz;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }
    
    if(phiy != 0) {
        /* rotate around y axis */
        rxx= cos(phiy); rxy= 0;         rxz= sin(phiy);
        //ryx= 0;         ryy= 1;         ryz= 0;
        rzx=-sin(phiy); rzy= 0;         rzz= cos(phiy);
        
        rotated_x = new_x*rxx + new_y*rxy + new_z*rxz;
        rotated_y = new_y;
        rotated_z = new_x*rzx + new_y*rzy + new_z*rzz;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }
    
    if(phiz != 0){
        /* rotate around z axis */
        rxx= cos(phiz); rxy=-sin(phiz); rxz= 0;
        ryx= sin(phiz); ryy= cos(phiz); ryz= 0;
        //rzx= 0;         rzy= 0;         rzz= 1;
        
        rotated_x = new_x*rxx + new_y*rxy ;
        rotated_y = new_x*ryx + new_y*ryy;
        rotated_z = new_z;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }
    
    new[1]=new_x;
    new[2]=new_y;
    new[3]=new_z;
    
    return new;
}



/* rotate a point about a unit vector axis */
double *rotate_axis(double *v, double *new, double *axis, double phi) {

    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double dot = (axis[1]*v[1]+axis[2]*v[2]+axis[3]*v[3])*(1.0-cosphi);

    new[1] = axis[1]*dot+v[1]*cosphi+(-axis[3]*v[2]+axis[2]*v[3])*sinphi;
    new[2] = axis[2]*dot+v[2]*cosphi+(+axis[3]*v[1]-axis[1]*v[3])*sinphi;
    new[3] = axis[3]*dot+v[3]*cosphi+(-axis[2]*v[1]+axis[1]*v[2])*sinphi;

    return new;
}




/* returns a unit vector in a random direction in arguments dx,dy,dz */
/* also returns a random magnitude within the unit sphere as a return value */
float uniform3Ddev(float *dx, float *dy, float *dz, long *seed)
{
    float ran1(long *idum);
    float dr;
    
    /* pick a random direction by cutting a sphere out of a cube */
    dr = 0;
    while(dr>1 || dr < 1e-2)
    {
        *dx = 2.1*(ran1(seed)-0.5);
        *dy = 2.1*(ran1(seed)-0.5);
        *dz = 2.1*(ran1(seed)-0.5);
        dr = sqrt(*dx**dx+*dy**dy+*dz**dz);
    }
    /* turn this into a unit vector */
    *dx/=dr;
    *dy/=dr;
    *dz/=dr;
    
    /* dx,dy,dz should now be a random unit vector */
    
    return dr;
}


float poidev(float xm, long *idum)
{
    float gammln(float xx);
    float ran1(long *idum);
    /* oldm is a flag for whether xm has changed since last call */
    static float sq,alxm,g,oldm=(-1.0);
    float em,t,y;
        
    /* routine below locks up for > 1e6 photons? */
    if (xm > 1.0e6) {
        return xm+sqrt(xm)*gaussdev(idum);
    }

    if (xm < 12.0) {
        /* use direct method: simulate exponential delays between events */
        if(xm != oldm) {
            /* xm is new, compute the exponential */
            oldm=xm;
            g=exp(-xm);
        }
        /* adding exponential deviates is equivalent to multiplying uniform deviates */
        /* final comparison is to the pre-computed exponential */
        em = -1;
        t = 1.0;
        do {
            ++em;
            t *= ran1(idum);
        } while (t > g);
    } else {
        /* Use rejection method */
        if(xm != oldm) {
            /* xm has changed, pre-compute a few things... */
            oldm=xm;
            sq=sqrt(2.0*xm);
            alxm=log(xm);
            g=xm*alxm-gammln(xm+1.0);
        }
        do {
            do {
                /* y is a deviate from a lorentzian comparison function */
                y=tan(M_PI*ran1(idum));
                /* shift and scale */
                em=sq*y+xm;
            } while (em < 0.0);		/* there are no negative Poisson deviates */
            /* round off to nearest integer */
            em=floor(em);
            /* ratio of Poisson distribution to comparison function */
            /* scale it back by 0.9 to make sure t is never > 1.0 */
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
        } while (ran1(idum) > t);
    }
        
    return em;
}


/* return gaussian deviate with rms=1 and FWHM = 2/sqrt(log(2)) */
float gaussdev(long *idum)
{
    float ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
        
    if (iset == 0) {
        /* no extra deviats handy ... */
        
        /* so pick two uniform deviates on [-1:1] */
        do {
            v1=2.0*ran1(idum)-1.0;
            v2=2.0*ran1(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0);
        /* restrained to the unit circle */
        
        /* apply Box-Muller transformation to convert to a normal deviate */
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;		/* we now have a spare deviate */
        return v2*fac;
    } else {
        /* there is an extra deviate in gset */
        iset=0;
        return gset;
    }
}


/* generate Lorentzian deviate with FWHM = 2 */
float lorentzdev(long *seed) {
    float ran1(long *idum);
 
    return tan(M_PI*(ran1(seed)-0.5));
}

/* return triangular deviate with FWHM = 1 */
float triangledev(long *seed) {
    float ran1(long *idum);
    float value;

    value = ran1(seed);
    if(value > 0.5){
        value = sqrt(2*(value-0.5))-1;
    }else{
        value = 1-sqrt(2*value);
    }

    return value;
}



float expdev(long *idum)
{
    float dum;
    
    do
    dum=ran1(idum);
    while( dum == 0.0);
    return -log(dum);
}



/* ln of the gamma function */
float gammln(float xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser = 1.000000000190015;
    for(j=0;j<=5;++j) ser += cof[j]/++y;
    
    return -tmp+log(2.5066282746310005*ser/x);
}





/* returns a uniform random deviate between 0 and 1 */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0 || !iy) {
        /* first time around.  don't want idum=0 */
        if(-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        
        /* load the shuffle table */
        for(j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if(*idum < 0) *idum += IM;
            if(j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    /* always start here after initializing */
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}


void polint(double *xa, double *ya, double x, double *y)
{
	double x0,x1,x2,x3;
	x0 = (x-xa[1])*(x-xa[2])*(x-xa[3])*ya[0]/((xa[0]-xa[1])*(xa[0]-xa[2])*(xa[0]-xa[3])); 
        x1 = (x-xa[0])*(x-xa[2])*(x-xa[3])*ya[1]/((xa[1]-xa[0])*(xa[1]-xa[2])*(xa[1]-xa[3]));
	x2 = (x-xa[0])*(x-xa[1])*(x-xa[3])*ya[2]/((xa[2]-xa[0])*(xa[2]-xa[1])*(xa[2]-xa[3]));
	x3 = (x-xa[0])*(x-xa[1])*(x-xa[2])*ya[3]/((xa[3]-xa[0])*(xa[3]-xa[1])*(xa[3]-xa[2]));
	*y = x0+x1+x2+x3;
}









/* read in multi-column text file to list of double arrays */
/* provide address of undeclared arrays on command line */
size_t read_text_file(char *filename, size_t nargs, ... )
{
    /* maximum of 10240-character lines? */
    char text[10240];
    char *token;
    const char delimiters[] = " \t,;:!";
    const char numberstuf[] = "0123456789-+.EGeg";

    unsigned long line,lines;
    unsigned long i,j;
    double value;
    double *data;
    double **pointer;
    va_list arglist;    
    FILE *infile = NULL;
    
    infile = fopen(filename,"r");
    if(infile == NULL) {
	perror("fopen()");
	return 0;
    }
    lines=0;
    while ( fgets ( text, sizeof text, infile ) != NULL ) {
	token = text;
	token += strspn(token,delimiters);
	if(strcmp(token,"\n")==0) {
	    //printf("blank\n");
	    continue;
	}
	++lines;
    }
    rewind(infile);

    /* allocate memory for arrays */
    va_start( arglist, nargs);
    for(i=0;i<nargs;++i){
	/* allocate the array */
	data = malloc((lines+10)*sizeof(double));
	/* initialize with missing number flags */
	for(j=0;j<lines+10;++j) {
	    data[j] = NAN;
	}
	/* get argument (pointer to pointer) */
	pointer = va_arg(arglist, double **);
	/* change the value of what the arg points to */
	*pointer = data;
	/* now the pointer provided as an argument points to
	something */
    }
    va_end(arglist);
        
    line = 0;
    while ( fgets ( text, sizeof text, infile ) != NULL ) { /* read a line */

	token = text;
	token += strspn(token,delimiters);
	if(strcmp(token,"\n")==0) {
	    //printf("blank\n");
	    continue;
	}
	i=0;
        va_start( arglist, nargs);
	do
        {
	    value=atof(token);
	    /* get argument */
	    pointer = va_arg(arglist, double **);
	    /* retrieve data array's address */
	    data = *pointer;
	    data[line] = value;

	    token += strspn(token,numberstuf);
	    if (strcmp(token,"\n")==0) continue;
	    token += strcspn(token,delimiters);
	    token += strspn(token,delimiters);
	    if (strcmp(token,"\n")==0) continue;

	    ++i;
	    if(i>=nargs) {
	        break;
	    }
	}
	while (strcmp(token,"\n")!=0) ;
	va_end(arglist);
 
//	printf("initializing:");
//        va_start( arglist, nargs);
//        for(i=0;i<nargs;++i){
//	    pointer = va_arg(arglist, double **);
//	    data = *pointer;
//	    printf(" %g",data[line]);
//        }
//        va_end(arglist);
//	printf("\n");

	++line;
    }
    fclose(infile);

    return lines;
}



/* measure magnitude of provided vector */
double magnitude(double *vector) {

    /* measure the magnitude */
    vector[0] = sqrt(vector[1]*vector[1]+vector[2]*vector[2]+vector[3]*vector[3]);

    return vector[0];
}

/* make provided vector a unit vector */
double unitize(double *vector, double *new_unit_vector) {
    double mag;

    /* measure the magnitude */
    mag = magnitude(vector);

    if(mag != 0.0){
    	/* normalize it */
	new_unit_vector[1]=vector[1]/mag;
	new_unit_vector[2]=vector[2]/mag;
	new_unit_vector[3]=vector[3]/mag;
    }
    else
    {
    	/* can't normalize, report zero vector */
    	new_unit_vector[0] = 0.0;
    	new_unit_vector[1] = 0.0;
    	new_unit_vector[2] = 0.0;
    	new_unit_vector[3] = 0.0;
    }
    return mag;
}

/* scale magnitude of provided vector */
double vector_scale(double *vector, double *new_vector, double scale) {

    new_vector[1] = scale*vector[1];
    new_vector[2] = scale*vector[2];
    new_vector[3] = scale*vector[3];
    
    return magnitude(new_vector);
}

/* enforce magnitude of provided vector */
double vector_rescale(double *vector, double *new_vector, double new_magnitude) {
    double oldmag;

    oldmag = magnitude(vector);
    if(oldmag <= 0.0) oldmag = 1.0;
    new_vector[1] = new_magnitude/oldmag*vector[1];
    new_vector[2] = new_magnitude/oldmag*vector[2];
    new_vector[3] = new_magnitude/oldmag*vector[3];

    return magnitude(new_vector);
}

/* difference between two given vectors */
double vector_diff(double *vector, double *origin_vector, double *new_vector) {

    new_vector[1] = vector[1]-origin_vector[1];
    new_vector[2] = vector[2]-origin_vector[2];
    new_vector[3] = vector[3]-origin_vector[3];
    return magnitude(new_vector);
}


/* vector cross product where vector magnitude is 0th element */
double *cross_product(double *x, double *y, double *z) {
    z[1] = x[2]*y[3] - x[3]*y[2];
    z[2] = x[3]*y[1] - x[1]*y[3];
    z[3] = x[1]*y[2] - x[2]*y[1];
    z[0] = 0.0;

    return z;
}
/* vector inner product where vector magnitude is 0th element */
double dot_product(double *x, double *y) {
    return x[1]*y[1]+x[2]*y[2]+x[3]*y[3];
}


/* polarization factor */
double polarization_factor(double kahn_factor, double *incident, double *diffracted, double *axis)
{
    double cos2theta,cos2theta_sqr,sin2theta_sqr;
    double twotheta,psi;
    double E_in[4];
    double B_in[4];
    double E_out[4];
    double B_out[4];
    
    unitize(incident,incident);
    unitize(diffracted,diffracted);
    unitize(axis,axis);
    
    /* component of diffracted unit vector along incident beam unit vector */
    cos2theta = dot_product(incident,diffracted);
    cos2theta_sqr = cos2theta*cos2theta;
    sin2theta_sqr = 1-cos2theta_sqr;

    if(kahn_factor != 0.0){
        /* tricky bit here is deciding which direciton the E-vector lies in for each source
           here we assume it is closest to the "axis" defined above */

        /* cross product to get "vertical" axis that is orthogonal to the cannonical "polarization" */
	cross_product(axis,incident,B_in);
        /* make it a unit vector */
        unitize(B_in,B_in);

        /* cross product with incident beam to get E-vector direction */
	cross_product(incident,B_in,E_in);
        /* make it a unit vector */
        unitize(E_in,E_in);

        /* get components of diffracted ray projected onto the E-B plane */
        E_out[0] = dot_product(diffracted,E_in);
        B_out[0] = dot_product(diffracted,B_in);

        /* compute the angle of the diffracted ray projected onto the incident E-B plane */
        psi = -atan2(B_out[0],E_out[0]);
    }

    /* correction for polarized incident beam */
    return 0.5*(1.0 + cos2theta_sqr - kahn_factor*cos(2*psi)*sin2theta_sqr);
}




SMVinfo GetFrame(char *filename)
{
    char *string;
    SMVinfo frame;
    char *byte_order;
    unsigned short int tempint;

    typedef union
    {
	unsigned char string[2];
	unsigned short integer;
    } TWOBYTES;
    TWOBYTES twobytes;
    twobytes.integer = 24954;
    

    /* determine byte order on this machine */
    if(0==strncmp((const char *) twobytes.string, "az", 2))
    {
	byte_order = "big_endian";
    }
    else
    {
	byte_order = "little_endian";
    }

    /* try to open the file... */
    frame.handle = fopen(filename, "rb");
    if(frame.handle != NULL)
    {
        /* just assume header will be 512 bytes?... */
        frame.header = calloc(1024,sizeof(char));
	if(! fread(frame.header, 512, 1, frame.handle))
	{
	    perror("SMV file header");
	    exit(9);
	}
	string = frame.header + 512;
        *string = (char) 0;

	/* remember the file name */
	frame.filename = calloc(strlen(filename)+10,sizeof(char));
	strcpy(frame.filename,filename);

	/* What kind of file is this? */
	if(0!=strncmp(frame.header, "{\nHEADER_BYTES=  512;\nDIM=2;\nBYTE_ORDER=", 12))
	{
	    /* probably not an ADSC frame */

	    /* inform the user */
	    printf("ERROR: %s does not look like an ADSC frame!\n", filename);
	    /* skip this file */
	    fclose(frame.handle);
	    
	    frame.handle = NULL;
	}
	else
	{
	    /* store the full header */
	    frame.header_size = (int) ValueOf("HEADER_BYTES",frame);
	    if(frame.header_size != 512)
	    {
		free(frame.header);
		fseek(frame.handle,0,SEEK_SET);
		frame.header = calloc(2*frame.header_size,sizeof(char));
		if(! fread(frame.header, frame.header_size, 1, frame.handle))
		{
		    perror("SMV file fread");
		    exit(9);
		}
		string = frame.header + frame.header_size;
	        *string = (char) 0;		
	    }

	    /* see if we will need to swap bytes */
	    string = (char *) strstr(frame.header, "BYTE_ORDER=")+11;
	    /* find last instance of keyword in the header */
	    while ((char *) strstr(string, "BYTE_ORDER=") != NULL)
	    {
		string = (char *) strstr(string, "BYTE_ORDER=")+11;
	    }
	    if(0==strncmp(byte_order, string, 10))
	    {
		frame.swap_bytes = FALSE;
	    }
	    else
	    {
		frame.swap_bytes = TRUE;
	    }

	    /* store a couple of things */
	    frame.width  = (int) ValueOf("SIZE1",frame);
	    frame.height = (int) ValueOf("SIZE2",frame);

	    if(frame.width == 0)
	    {
		/* try other formats? */
		frame.width = frame.height = (int) ValueOf("DETECTOR_DIMENSIONS",frame);
	    }

//	    frame.mmapdata = mmap(NULL,2*frame.width*frame.height+frame.header_size,PROT_READ,MAP_SHARED,fileno(frame.handle),0);
	    frame.mmapdata = calloc(2,frame.width*frame.height+frame.header_size);
	    if(frame.mmapdata == NULL)
	    {
		perror("calloc:");
	    }
	    fseek(frame.handle,0,SEEK_SET);
	    printf("reading %s\n",frame.filename);
	    if(! fread(frame.mmapdata,1,2*frame.width*frame.height+frame.header_size,frame.handle))
	    {
	        perror("SMV file fread");
	        exit(9);
	    }

	    printf("mmap(%s) = %p\n",frame.filename,frame.mmapdata);


	}
    }
    else
    {
	/* fopen() failed */
	perror("nonBragg");
    }
    
    return frame;
}

/* read floating-point values from keywords in an SMV header */
double ValueOf(const char *keyword, SMVinfo frame)
{
    double value;
    char *string;
    int keylen = strlen(keyword);

    /* start at the beginning */
    string = frame.header;

    /* find first instance of keyword in the header */
//    string = (char *) strstr(frame.header, keyword);
//    string = string + keylen;
    /* find last instance of keyword in the header */
    while ((char *) strstr(string, keyword) != NULL)
    {
	string = (char *) strstr(string, keyword)+keylen;
    }
    if(string == frame.header) return NAN;

    /* advance to just after the "=" sign */
    string = (char *) strstr(string, "=");
    if(string == NULL) return 0.0;
    ++string;

    value = atof(string);

    return value;
}


unsigned char *read_pgm5_bytes(char *filename,unsigned int *returned_width,unsigned int *returned_height)
{
    unsigned char test[512];
    unsigned char *array = NULL;
    FILE *handle = NULL;
    unsigned int width=0,height=0,maxvalue=0;

    handle = fopen(filename,"rb");
    if(handle)
    {
        if(! fread(test,512,1,handle))
	{
	    perror("PGM fread header");
	    exit(9);
	}
        if(strstr((const char *) test,"P5"))
        {
            /* PGM header: "P5<whitespace>width<whitespace>height<whitespace>maxvalue<single whitespace character>" */
            fseek(handle,3,SEEK_SET);
            if(! fscanf(handle," %u %u %u",&width,&height,&maxvalue))
	    {
	        perror("PGM fscanf");
	        exit(9);
	    }
            /* skip final single whitespsce character (first pixel could have value of "20") */
            fseek(handle,1,SEEK_CUR);
            array = calloc(sizeof(unsigned char),width*height);
            if(! fread(array,width,height,handle))
	    {
	        perror("PGM fread");
	        exit(9);
	    }
        }
        fclose(handle);
    }
    else
    {
        perror("PGM fopen");
    }

    *returned_width = width;
    *returned_height = height;
    return array;
}



int compare_float(const void *ptr1,const void *ptr2){
    int result = 0;
    float first,second;

    first = *( (float *) ptr1);
    second = *( (float *) ptr2);

    if(first < second) result = -1;
    if(first == second) result = 0;
    if(first > second) result =  1;

    return result;
}



#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
float fmedian(unsigned int n, float arr[])
{
    unsigned int i,j,k,l,ir,mid;
    float a,temp;

    l=1;
    ir=n;
    k=(n+1)/2;
//printf("n=%d; k=%d\n",n,k);

//for(i=1;i<=n;++i) printf("arr[%d]=%f\n",i,arr[i]);

    for(;;)
    {
	if(ir <= l+1)
	{
	    if(ir == l+1 && arr[ir] < arr[l])
	    {
		SWAP(arr[l],arr[ir]);
	    }
//for(i=1;i<=n;++i) printf("arr[%d]=%f\n",i,arr[i]);
	    return arr[k];
	} else {
	    mid=(l+ir) >> 1;
	    SWAP(arr[mid],arr[l+1]);
	    if(arr[l+1] > arr[ir])
	    {
		SWAP(arr[l+1],arr[ir]);
	    }
	    if(arr[l] > arr[ir])
	    {
		SWAP(arr[l],arr[ir]);
	    }
	    if(arr[l+1] > arr[l])
	    {
		SWAP(arr[l+1],arr[l]);
	    }
	    i=l+1;	// initialize pointers for partitioning
	    j=ir;	
	    a=arr[l];	// partitioning element
	    for(;;)	// innermost loop
	    {
		do i++; while(arr[i]<a);	// scan up to find element > a
		do j--; while(arr[j]>a);	// scan down to find element < a
		if( j < i ) break;		// pointers crossed, median is in between
		SWAP(arr[i],arr[j]);
	    }
	    arr[l]=arr[j];			// insert partitioning element
	    arr[j]=a;
	    if( j >= k ) ir=j-1;		// Keep partition that contains the median active
	    if( j <= k ) l=i;
	}
    }
}


float fmedian_with_rejection(unsigned int n, float arr[],float sigma_cutoff, float *final_mad, int *final_n)
{
    float median_value;
    int i,orig_n,reject,worst,done;
    float min_frac,sum,deviate,mad,worst_deviate,temp;

    orig_n = n;
    min_frac = 0.7;

    done = 0;
    while(! done)
    {
	/* compute the median (centroid) value */
	median_value = fmedian(n,arr);

	/* now figure out what the mean absolute deviation from this value is */
	mad = fmedian_absolute_deviation(n,arr,median_value);
	//if(flag) printf("mad = %f\n",mad);

	done = 1;
	/* reject all outliers */
	for(i=1;i<=n;++i)
	{
	    /* reject positive and negative outliers */
	    deviate = fabs(arr[i]-median_value);
	    if(deviate > sigma_cutoff*mad)
	    {
	        /* needs to go */
	        /* move value at the end of the array to this "reject" and then shorten the array */
	        //if(flag) printf("rejecting arr[%d] = %f (%f)\n",i,arr[i],deviate);
	        //arr[worst]+=10000;
	        if(i != n)
	        {
		    //temp=arr[worst];
		    arr[i] = arr[n];
		    //arr[n]=temp;
		}
		--n;
		done = 0;
	    }
	}
    }

    /* basically three return values */
    *final_mad = mad;
    *final_n = n;
    return median_value;
}

/* note: there must be 2*n elements in this array! */
float fmedian_absolute_deviation(unsigned int n, float arr[], float median_value)
{
    int i;
    for(i=1;i<=n;++i)
    {
	arr[i+n] = fabs(arr[i]-median_value);
    }

    return fmedian(n,arr+n);
}




/* this function keeps track of outliers by swapping them to the end of the array */
/* counting starts at 0 and "points" is the number of points */
float fmean_with_rejection(unsigned int starting_points, float arr[], float sigma_cutoff, float *final_rmsd, int *final_n)
{
    int points,n,i;
    int rejection,worst;
    float temp,sum,avg,sumd,rmsd,deviate,worst_deviate;

    points=starting_points;
    rejection = 1;
    while ( rejection && points>starting_points/2.0 )
    {
        /* find the mean and rms deivation */
        sum = sumd = 0.0;
        for(i=0;i<points;++i)
        {
	    sum+=arr[i];
        }
        avg=sum/points;
	worst=-1;
	worst_deviate=0.0;
        for(i=0;i<points;++i)
        {
	    deviate=fabs(arr[i]-avg);
	    if(deviate > worst_deviate)
	    {
		worst=i;
		worst_deviate=deviate;
	    }
	    sumd+=deviate*deviate;
        }
        rmsd=sqrt(sumd/points);

	rejection=0;
	if(worst_deviate>sigma_cutoff*rmsd)
	{
	    /* we have a reject! */
	    rejection=1;

	    /* move it to end of the array and forget about it */
	    SWAP(arr[worst],arr[points]);
	    --points;
	}
    }

    *final_rmsd = rmsd;
    *final_n = points;
    return avg;
}


