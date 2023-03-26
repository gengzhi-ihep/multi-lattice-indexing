/* convert ideal pixel intensities into noisy pixels                                            -James Holton           6-9-17

example:

gcc -O -o noisify noisify.c -lm

./noisify -bin floatimage.bin -distance 100 -detsize 100 -pixel 0.1 \
  -scale 1 -readout 3 -flicker 0.02 -calibration 0.03

wavelength (lambda) should be provided in Angstrom
detector distance, detsize and pixel size in mm
the -scale value is multiplied by every value found in floatimage.bin before use

floatimage.bin should be a binary "dumpfile" consisting of the proper number of 4-byte
"float" numbers on the current architecture.  These numbers should be in "photons/pixel" scale.
The nearBragg and fastBragg programs can be used to generate it.

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

typedef enum { UNKNOWN, FIBER, GAUSS
 } psf_type;
float *apply_psf(float *inimage, int xpixels, int ypixels, psf_type psftype, double fwhm_pixels, int psf_radius);

/* analytic integral of a Gaussian */
double ngauss2D(double x, double y, double fwhm);
double ngauss2D_integ(double x, double y);
double ngauss2D_pixel(double x,double y,double pix);
double integrate_gauss_over_pixel(double x, double y, double fwhm, double pix);

/* analytic integral of fiber PSF function */
double fiber2D_integ(double x, double y,double g);
double fiber2D_pixel(double x,double y,double g, double pix);
double integrate_fiber_over_pixel(double x, double y, double g, double pix);



char *floatfilename = "floatimage.bin\0";
FILE *floatfile = NULL;
char *headerfilename = NULL;
SMVinfo headerfile;
char *intfilename = "intimage.img\0";
char *pgmfilename = "image.pgm\0";
char *noisefilename = "noiseimage.img\0";
FILE *outfile = NULL;

int main(int argc, char** argv)
{

    /* detector parameters used to make the header */
    /* assumed to be the same as those used to call nearBragg/fastBragg!  */

    double distance = 100.0e-3;
    double detsize_x = 102.4e-3;
    double detsize_y = 102.4e-3;
    double pixel = 0.1e-3;
    double Xdet,Ydet,Xbeam=-1e99,Ybeam=-1e99,Rdet;
    int xpixel,ypixel,xpixels=0,ypixels=0,pixels;
    double lambda = 1;

    psf_type psftype = UNKNOWN;
    float psf_fwhm = 46e-6;
    int psf_radius = 0;
    int x0,y0,x,y,dx,dy;
    float rsq,temp;

    int n,i,j;
    float *floatimage,*photonimage,*psfimage,*spare;
    unsigned short int *int16image;
    unsigned int *int32image;
    unsigned char *pgmimage;

    double test,sum,photons,photons0,adu;
    double readout_noise=0.0, flicker_noise=0.0;
    double calibration_noise=0.03;
    double adc_offset = 40.0;
    double quantum_gain = 1.0;
    int overloads = 0;

    int calculate_noise = 1;
    int write_pgm = 1;

    double phi0 = 0, osc = 1;

    /* Thomson cross section */
    double r_e_sqr = 7.94079248018965e-30;
    /* incident x-ray fluence in photons/m^2   default equivalent to unity
        that is, one electron will scatter 1 ph/SR after a fluence of 1.26e29 ph/m^2
        this places the input file on a photons/pixel scale */
    double fluence = 125932015286227086360700780544.0;
    /* arbitrary "photon scale" applied before calculating noise, default is unity */
    double photon_scale = 1.0;
    double intfile_scale;

    double I;
    double max_I = 0.0;

    long seed;
    long calib_seed = 123456789;

    seed = -time((time_t *)0);
//    printf("GOTHERE seed = %u\n",seed);


    /* check argument list */
    for(i=1; i<argc; ++i)
    {
        if(argv[i][0] == '-')
        {
            /* option specified */
            if(strstr(argv[i], "-lambda") && (argc > (i+1)))
            {
                /* copy directly into image header */
                lambda = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-Xbeam") && (argc > (i+1)))
            {
                Xbeam = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-Ybeam") && (argc > (i+1)))
            {
                Ybeam = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-distance") && (argc > (i+1)))
            {
                distance = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detsize") && (strlen(argv[i]) == 8) && (argc > (i+1)))
            {
                detsize_x = atof(argv[i+1])/1000.0;
                detsize_y = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detsize_x") && (argc > (i+1)))
            {
                detsize_x = atof(argv[i+1])/1000.0;
            }
             if(strstr(argv[i], "-detsize_y") && (argc > (i+1)))
            {
                detsize_y = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detpixels") && (strlen(argv[i]) == 10) && (argc > (i+1)))
            {
                xpixels = ypixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detpixels_x") && (argc > (i+1)))
            {
                xpixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detpixels_y") && (argc > (i+1)))
            {
                ypixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-pixel") && (argc > (i+1)))
            {
                pixel = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-psf") && (strlen(argv[i]) == 4) && (argc >= (i+1)))
            {
                psftype = UNKNOWN;
                if(strstr(argv[i+1],"gauss")) psftype = GAUSS;
                if(strstr(argv[i+1],"fiber")) psftype = FIBER;
                if(psftype == UNKNOWN) printf("WARNING: unknown psf type: %s\n",argv[i+1]);
            }
            if(strstr(argv[i], "-psf_rad") && (argc > (i+1)))
            {
                psf_radius = atof(argv[i+1]);
            }
            if((strstr(argv[i], "-psf_si") || strstr(argv[i], "-psf_fw") || strstr(argv[i], "-psf_wi")) && (argc > (i+1)))
            {
                psf_fwhm = atof(argv[i+1])/1e6;
            }
            if(strstr(argv[i], "-fluence") && (argc > (i+1)))
            {
                fluence = atof(argv[i+1]);
            }
            if((strstr(argv[i], "-floatfile") || strstr(argv[i], "-floatimage") || strstr(argv[i], "-bin")) && (argc > (i+1)))
            {
                floatfilename = argv[i+1];
                floatfile = fopen(floatfilename,"r");
            }
            if(strstr(argv[i], "-header") && (argc > (i+1)))
            {
                headerfilename = argv[i+1];
            }
            if((strstr(argv[i], "-pgmfile") || strstr(argv[i], "-pgmimage")) && (argc > (i+1)))
            {
                pgmfilename = argv[i+1];
            }
            if((strstr(argv[i], "-intfile") || strstr(argv[i], "-intimage")) && (argc > (i+1)))
            {
                intfilename = argv[i+1];
            }
            if((strstr(argv[i], "-noisefile") || strstr(argv[i], "-noiseimage")) && (argc > (i+1)))
            {
                noisefilename = argv[i+1];
            }
            if(strstr(argv[i], "-nonoise") )
            {
                /* turn off noise */
                calculate_noise = 0;
            }
            if(strstr(argv[i], "-nopgm") )
            {
                write_pgm = 0;
            }
            if((strstr(argv[i], "-readout") || strstr(argv[i], "-readnoi")) && (argc > (i+1)))
            {
                readout_noise = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-flicker") && (argc > (i+1)))
            {
                flicker_noise = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-calibration") && (argc > (i+1)))
            {
                calibration_noise = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-scale") && (argc > (i+1)))
            {
                photon_scale = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-seed") && (argc > (i+1)))
            {
                seed = -atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-calib_seed") && (argc > (i+1)))
            {
                calib_seed = -atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-adc") && (argc > (i+1)))
            {
                adc_offset = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-gain") && (argc > (i+1)))
            {
                quantum_gain = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-phi") && (argc > (i+1)))
            {
                phi0 = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-osc") && (argc > (i+1)))
            {
                osc = atof(argv[i+1]);
            }
        }
    }

    printf("noisify - add noise to pixels - James Holton 2-16-16\n");

    if(floatfile == NULL){
        printf("usage: noisify -floatfile floatimage.bin\n");
        printf("options:\n");\
        printf("\tfloatimage.bin\t nearBragg-style binary dump file\n");
        printf("\t-scale\tscale factor to put floatimage.bin in photons/pixel\n");
        printf("\t-gain\tpixel units per photon\n");
        printf("\t-readout_noise\tgaussian noise added to every pixel\n");
        printf("\t-flicker\t fractional 1/f noise in source\n");
        printf("\t-calibration\t static fractional error per pixel\n");
        printf("\t-calib_seed\t change seed for calibration error\n");
        printf("\t-seed\t specify seed for all non-calibration errors\n");
        printf("\t-gain\t pixel units per photon\n");
        printf("\t-adc\t offset added to each pixel after noise\n");
        printf("\t-distance\t distance from origin to detector center in mm\n");
        printf("\t-detsize\t detector size in mm\n");
        printf("\t-pixel\t detector pixel size in mm\n");
        printf("\t-psf gauss|fiber\t point spread function type (gaussian or fiber)\n");
        printf("\t-psf_fwhm\t point spread function size in um\n");
        printf("\t-psf_radius\t radius to render PSF in pixels (default automatic)\n");
        printf("\t-lambda\t incident x-ray wavelength in Angstrom\n");
        printf("\t-intfile\t name of smv-formatted output file (arbitrary scale)\n");
        printf("\t-pgmfile\tname of pgm-formatted output file (arbitrary scale)\n");
        printf("\t-noisefile\t name of smv-formatted output file (with noise)\n");
        printf("\t-Xbeam\t image X coordinate of direct-beam spot (mm)\n");
        printf("\t-Ybeam\t image Y coordinate of direct-beam spot (mm)\n");
        printf("\t-header\t import 512-byte header from specified SMV file\n");
exit(9);
    }

    /* count how much data we got */
    fseek(floatfile,0,SEEK_END);
    n = ftell(floatfile);
    rewind(floatfile);
    pixels = n/sizeof(float);

    if(headerfilename != NULL)
    {
        printf("taking header from %s\n",headerfilename);
        /* frame handling routines */
        headerfile = GetFrame(headerfilename);
        if(headerfile.header_size > 0) {
            xpixels = headerfile.width;
            ypixels = headerfile.height;
            pixels = xpixels*ypixels;
            test = ValueOf("PIXEL_SIZE",headerfile);
            if(! isnan(test)) pixel = test/1000.0;
            detsize_x = pixel*xpixels;
            detsize_y = pixel*ypixels;
            test = ValueOf("DISTANCE",headerfile);
            if(! isnan(test)) distance = test/1000.0;
//          test = ValueOf("CLOSE_DISTANCE",headerfile);
//          if(! isnan(test)) close_distance = test/1000.0;
            test = ValueOf("WAVELENGTH",headerfile);
            if(! isnan(test)) lambda = test/1e10;
            test = ValueOf("BEAM_CENTER_X",headerfile);
            if(! isnan(test)) Xbeam = test/1000.0;
            test = ValueOf("BEAM_CENTER_Y",headerfile);
            if(! isnan(test)) Ybeam = detsize_y - test/1000.0;
//          test = ValueOf("ORGX",headerfile);
//          if(! isnan(test)) ORGX = test;
//          test = ValueOf("ORGY",headerfile);
//          if(! isnan(test)) ORGY = test;
//          test = ValueOf("PHI",headerfile);
//          if(! isnan(test)) phi0 = test/RTD;
//          test = ValueOf("OSC_RANGE",headerfile);
//          if(! isnan(test)) osc = test/RTD;
//          test = ValueOf("TWOTHETA",headerfile);
//          if(! isnan(test)) twotheta = test/RTD;
        }
    }

    /* other sensibe defaults */
    if(! xpixels && ! ypixels) {
        /* hmm... guess? */
        printf("WARNING: guessing xy pixel dimensions.\n");
        xpixels = sqrt(pixels);
        ypixels = pixels/xpixels;
        while( pixels != xpixels*ypixels && xpixels > 0 )
        {
            --xpixels;
            ypixels = pixels/xpixels;
        }
        if( pixels != xpixels*ypixels) {
             xpixels = pixels;
             ypixels = 1;
        }
    }
    if(xpixels && ! ypixels) {
        ypixels = pixels/xpixels;
    }
    if(! xpixels && ypixels) {
        xpixels = pixels/ypixels;
    }

    /* finalize detector size */
    if(xpixels) {
        detsize_x = pixel*xpixels;
    }
    else
    {
        xpixels = ceil(detsize_x/pixel-0.5);
    }
    if(ypixels) {
        detsize_y = pixel*ypixels;
    }
    else
    {
        ypixels = ceil(detsize_y/pixel-0.5);
    }
    pixels = xpixels*ypixels;

    /* allocate memory */
    floatimage = calloc(pixels+10,sizeof(float));
    photonimage = calloc(pixels+10,sizeof(float));
    int16image = calloc(pixels+10,sizeof(unsigned short int));
    int32image = calloc(pixels+10,sizeof(unsigned int));
    if(write_pgm) pgmimage   = calloc(pixels+10,sizeof(unsigned char));

    printf("importing %d pixel intensites: %s\n",pixels,floatfilename);
    if(! fread(floatimage,pixels,sizeof(float),floatfile))
    {
        perror("reading input file");
        exit(9);
    }
    fclose(floatfile);

    /* default to middle of detector unless specified earlier */
    if(Xbeam <= -1e99) Xbeam = detsize_x/2.0;
    if(Ybeam <= -1e99) Ybeam = detsize_y/2.0;

    if(calculate_noise == 0)
    {
        calibration_noise = 0;
        readout_noise = 0;
        flicker_noise = 0;
    }

    printf("  distance=%g detsize=%gx%g  pixel=%g meters (%dx%d pixels)\n",distance,detsize_x,detsize_y,pixel,xpixels,ypixels);
    printf("  Xbeam=%g Ybeam=%g\n",Xbeam,Ybeam);
    if(psftype == GAUSS) printf("  Gaussian PSF fwhm = %g um ",psf_fwhm*1e6);
    if(psftype == FIBER) printf("  fiber PSF fwhm = %g um ",psf_fwhm*1e6);
    if(psftype != UNKNOWN && psf_radius == 0) printf("  with automatic rendering radius\n");
    if(psftype != UNKNOWN && psf_radius >= 0) printf("  with rendering radius: %d\n",psf_radius);
    printf("  seed: %ld\n",seed);
    printf("  calibration noise seed: %ld\n",calib_seed);
    printf("  calibration_noise = %g %%\n",calibration_noise*100);
    printf("  input file scale = %g\n",photon_scale);
    printf("  readout_noise = %g ADU\n",readout_noise);
    printf("  flicker_noise = %g %%\n",flicker_noise*100);
    printf("  quantum_gain = %g ADU/photon\n",quantum_gain);
    printf("  adc_offset = %g ADU\n",adc_offset);


    printf("\n");


    /* put on photon scale first */
    max_I = 0.0;
    for(i=0;i<pixels;++i)
    {
        I = floatimage[i];
        if(max_I < I) max_I = I;
        if(I < 0.0) printf("WARNING: negative intensity in %s: %g\n",floatfilename,I);

        /* convert into photons/pixel (no change unless user specified fluence) */
        photonimage[i] = (fluence*r_e_sqr)*photon_scale*I;
    }
    printf("maximum value in input file: %g ( %g on photon scale)\n",max_I,max_I*photon_scale*fluence*r_e_sqr);


    /* do PSF on noiseless image only if it won't be available in the noise image */
    if(calculate_noise == 0 && psftype != UNKNOWN && psf_fwhm > 0.0)
    {
        /* run the blurring routine */
        printf("  applying PSF to noiseless image width = %g pixels\n",psf_fwhm/pixel);
        psfimage = apply_psf(photonimage, xpixels, ypixels, psftype, psf_fwhm/pixel, psf_radius);

        /* we won't be using photonimage data again. but what if apply_psf didn't calloc? */
//      free(photonimage);
        photonimage = psfimage;
    }


    /* output noiseless image as ints */
    for(i=0;i<pixels;++i)
    {
        /* convert noiseless photons/pixel into area detector units */
        adu = photonimage[i]*quantum_gain+adc_offset;
        if(adu > 65535.0) adu = 65535.0;
        int16image[i] = (unsigned short int) ( adu );
        //printf("%.50g %d\n",adu,int16image[i]);
    }
    printf("writing %s as %d %lu-byte integers\n",intfilename,pixels,sizeof(unsigned short int));
    outfile = fopen(intfilename,"wb");
    if(headerfilename != NULL)
    {
        /* use the provided header if possible */
        fwrite(headerfile.header,1,headerfile.header_size,outfile);
    }
    else
    {
        /* make up our own header */
        fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=little_endian;\nTYPE=unsigned_short;\n");
        fprintf(outfile,"SIZE1=%d;\nSIZE2=%d;\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",xpixels,ypixels,pixel*1000.0,distance*1000.0);
        fprintf(outfile,"WAVELENGTH=%g;\nBEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",lambda,Xbeam*1000.0,(detsize_y-Ybeam)*1000);
        fprintf(outfile,"PHI=%g;\nOSC_START=%g;\nOSC_RANGE=%g;\n",phi0,phi0,osc);
        fprintf(outfile,"TIME=%g;\n",osc);
        fprintf(outfile,"DETECTOR_SN=000;\n");
        fprintf(outfile,"BEAMLINE=fake;\n");
        fprintf(outfile,"}\f");
        while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    }
    fwrite(int16image,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);


    if(write_pgm)
    {
        /* output as pgm */
        for(i=0;i<pixels;++i){
            test = int16image[i];
            if(test > 255.0) test = 255.0;
            pgmimage[i] = (unsigned char) ( test );
            //printf("%d %d = %d\n",xpixel,ypixel,pgmimage[i]);
        }

        printf("writing %s as %lu-byte integers\n",pgmfilename,sizeof(unsigned char));
        outfile = fopen(pgmfilename,"wb");
        fprintf(outfile, "P5\n%d %d\n", xpixels, ypixels);
        fprintf(outfile, "# pixels scaled by %lg\n", 1.0);
        fprintf(outfile, "255\n");
        fwrite(pgmimage,sizeof(unsigned char),pixels,outfile);
        fclose(outfile);
    }

    /* quit now if there is nothing else to do */
    if(calculate_noise == 0){
        return 0;
    }

    /* simulate noise */
    sum = 0.0;
    for(i=0;i<pixels;++i){

        /* ideal photons/pixel */
        photons0 = photonimage[i];

        /* simulate 1/f noise in source */
        if(flicker_noise > 0.0){
            photons0 *= ( 1.0 + flicker_noise * gaussdev( &seed ) );
        }
        /* calibration is same from shot to shot, so use different seed */
        if(calibration_noise > 0.0){
            photons0 *= ( 1.0 + calibration_noise * gaussdev( &calib_seed ) );
        }
        /* simulate photon-counting error (assume calibration error is loss of photons, not electrons) */
        photonimage[i] = poidev( photons0, &seed );

        /* accumulate number of photons */
        sum += photonimage[i];
    }

    /* now that we have photon count at each point, implement any PSF */
    if(psftype != UNKNOWN && psf_fwhm > 0.0)
    {
        /* report on sum before the PSF is applied */
        printf("%.0f photons on noise image before PSF\n",sum);
        /* start with a clean slate */
        printf("  applying PSF width = %g um\n",psf_fwhm*1e6);
        psfimage = apply_psf(photonimage, xpixels, ypixels, psftype, psf_fwhm/pixel, psf_radius);

        /* from now on, this is the "photonimage", or singal that is subject to read noise */
//      free(photonimage);
        photonimage = psfimage;
    }


    sum = 0;
    overloads = 0;
    for(i=0;i<pixels;++i){
        sum += photonimage[i];

        /* convert photon signal to pixel units */
        adu = photonimage[i]*quantum_gain + adc_offset;

        /* readout noise is in pixel units? */
        if(readout_noise > 0.0){
            adu += readout_noise * gaussdev( &seed );
        }

        if(adu > 65535.0) {
            adu = 65535.0;
            ++overloads;
        }
        int16image[i] = (unsigned short int) adu;
//      printf("pixel %d = %d\n",i,int16image[i]);
    }
    printf("%.0f photons on noise image (%d overloads)\n",sum,overloads);

    printf("writing %s as %lu-byte integers\n",noisefilename,sizeof(unsigned short int));
    outfile = fopen(noisefilename,"wb");
    if(headerfilename != NULL)
    {
        /* use provided header if we have one */
        fwrite(headerfile.header,1,headerfile.header_size,outfile);
    }
    else
    {
        /* make up our own header */
        fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=little_endian;\nTYPE=unsigned_short;\n");
        fprintf(outfile,"SIZE1=%d;\nSIZE2=%d;\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",xpixels,ypixels,pixel*1000.0,distance*1000.0);
        fprintf(outfile,"WAVELENGTH=%g;\nBEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",lambda,Xbeam*1000.0,(detsize_y-Ybeam)*1000);
        fprintf(outfile,"PHI=%g;\nOSC_START=%g;\nOSC_RANGE=%g;\n",phi0,phi0,osc);
        fprintf(outfile,"TIME=%g;\n",osc);
        fprintf(outfile,"DETECTOR_SN=000;\n");
        fprintf(outfile,"BEAMLINE=fake;\n");
        fprintf(outfile,"}\f");
        while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    }
    fwrite(int16image,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);

    return 0;
}


/* 2D Gaussian integral=1 */
double ngauss2D(double x, double y, double fwhm)
{
    return log(16.)/M_PI*fwhm*fwhm*exp(-log(16.)*((x*x+y*y)/(fwhm*fwhm) ));
}

/* integral of Gaussian fwhm=1 integral=1 */
double ngauss2D_integ(double x, double y)
{
    return 0.125*(erf(2*x*sqrt(log(2.)))*erf(y*sqrt(log(16.)))*sqrt(log(16.)/log(2.)));
}

/* unit volume integrated over a pixel, fwhm = 1 */
double ngauss2D_pixel(double x,double y,double pix)
{
    return ngauss2D_integ(x+pix/2.,y+pix/2.)-ngauss2D_integ(x+pix/2.,y-pix/2.)-ngauss2D_integ(x-pix/2.,y+pix/2.)+ngauss2D_integ(x-pix/2.,y-pix/2.);
}

double integrate_gauss_over_pixel(double x, double y, double fwhm, double pix)
{
    return ngauss2D_pixel(x/fwhm,y/fwhm,pix/fwhm);
}


double fiber2D(double x,double y,double g)
{
    /* g/(2*pi)*(g**2+x**2+y**2)**(-3/2) */
    double temp;
    temp = sqrt(g*g+x*x+y*y);
    if(temp <= 0.0) return 0.0;
    return g/2.0/M_PI/temp/temp/temp;
}
double fiber2D_integ(double x,double y,double g)
{
    return atan((x*y)/(g*sqrt(g*g + x*x + y*y)))/2.0/M_PI;
}
double fiber2D_pixel(double x,double y,double g,double pix)
{
  return fiber2D_integ(x+pix/2.,y+pix/2.,g)-fiber2D_integ(x+pix/2.,y-pix/2.,g)-fiber2D_integ(x-pix/2.,y+pix/2.,g)+fiber2D_integ(x-pix/2.,y-pix/2.,g);
}
double integrate_fiber_over_pixel(double x, double y, double g, double pix)
{
    return fiber2D_pixel(x,y,g,pix);
}


/* function for applying the PSF, returns NEW image that is blurred version of input */
float *apply_psf(float *inimage, int xpixels, int ypixels, psf_type psftype, double fwhm_pixels, int user_psf_radius)
{
    double max_I;
    float *outimage;
    double *kernel;
    int x0,y0,x,y,dx,dy;
    double g,rsq;
    double photon_noise,lost_photons=0.0,total_lost_photons=0.0;
    int pixels,maxwidth,kernel_size,psf_radius;
    int i,j,k;
    double photonloss_factor = 10.0;

    /* convert fwhm to "g" distance : fwhm = sqrt((2**(2./3)-1))/2*g */
    g = fwhm_pixels * 0.652383013252053;

    if(psftype == UNKNOWN)
    {
        printf("ERROR: unknown PSF type\n");
        return inimage;
    }

    pixels = xpixels*ypixels;
    if(pixels == 0)
    {
        printf("ERROR: apply_psf image has zero size\n");
        return inimage;
    }

    if(fwhm_pixels <= 0.0)
    {
        printf("WARNING: apply_psf function has zero size\n");
        return inimage;
    }

    /* start with a clean slate */
    outimage = calloc(pixels+10,sizeof(float));

    psf_radius = user_psf_radius;
    if(psf_radius <= 0)
    {
        /* auto-select radius */

        /* preliminary stats */
        max_I = 0.0;
        for(i=0;i<pixels;++i)
        {
            /* optionally scale the input file */
            if(max_I < inimage[i]) max_I = inimage[i];
        }
        printf("  maximum input photon/pixel: %g\n",max_I);

        if(max_I<=0.0)
        {
            /* nothing to blur */
            printf("WARNING: no photons, PSF skipped\n");
            return outimage;
        }

        /* at what level will an error in intensity be lost? */
        photon_noise = sqrt(max_I);
        lost_photons = photon_noise/photonloss_factor;

        if(psftype == GAUSS)
        {
            /* calculate the radius beyond which only 0.5 photons will fall */
            psf_radius = 1+ceil( sqrt(-log(lost_photons/max_I)/log(4.)/2.)*fwhm_pixels );
            printf("  auto-selected psf_radius = %d pixels\n",psf_radius);
        }
        if(psftype == FIBER)
        {
            /* calculate the radius r beyond which only 0.5 photons will fall */
            /* r = sqrt((g*(max_I/0.5))**2-g**2)
                 ~ 2*g*max_I */
            psf_radius = 1+ceil( g*(max_I/lost_photons)  );
            printf("  auto-selected psf_radius = %d pixels\n",psf_radius);
        }
        if(psf_radius == 0) psf_radius = 1;
    }
    /* limit psf kernel to be no bigger than 4x the input image */
    maxwidth = xpixels;
    if(ypixels > maxwidth) maxwidth = ypixels;
    if(psf_radius > maxwidth) psf_radius = maxwidth;
    kernel_size = 2*psf_radius+1;

    /* now alocate enough space to store the PSF kernel image */
    kernel = calloc(kernel_size*kernel_size,sizeof(double));
    if(kernel == NULL)
    {
        perror("apply_psf: could not allocate memory for PSF kernel");
        exit(9);
    }

    /* cache the PSF in an array */
    for(dy=-psf_radius;dy<=psf_radius;++dy)
    {
        for(dx=-psf_radius;dx<=psf_radius;++dx)
        {
            rsq = dx*dx+dy*dy;
            if(rsq > psf_radius*psf_radius) continue;

            /* this could be more efficient */
            k = kernel_size*(kernel_size/2+dy)+kernel_size/2+dx;


            if( psftype == GAUSS ) {
                kernel[k] = integrate_gauss_over_pixel(dx,dy,fwhm_pixels,1.0);
            }
            if( psftype == FIBER ) {
                kernel[k] = integrate_fiber_over_pixel(dx,dy,g,1.0);
            }
        }
    }

    /* implement PSF  */
    for(i=0;i<pixels;++i)
    {
        x0 = i%xpixels;
        y0 = (i-x0)/xpixels;

        /* skip if there is nothing to add */
        if(inimage[i] <= 0.0) continue;

        if(user_psf_radius != 0)
        {
            psf_radius = user_psf_radius;
        }
        else
        {
            /* at what level will an error in intensity be lost? */
            photon_noise = sqrt(inimage[i]);
            lost_photons = photon_noise/photonloss_factor;

            if(psftype == GAUSS)
            {
                /* calculate the radius beyond which only 0.5 photons will fall
                   r = sqrt(-log(lost_photons/total_photons)/log(4)/2)*fwhm */
                psf_radius = 1+ceil( sqrt(-log(lost_photons/inimage[i])/log(16.))*fwhm_pixels );
//              printf("  auto-selected psf_radius = %d pixels\n",psf_radius);
            }
            if(psftype == FIBER)
            {
                /* calculate the radius beyond which only 0.5 photons will fall
                   r = sqrt((g*(total_photons/lost_photons))**2-g**2)
                     ~ g*total_photons/lost_photons */
                psf_radius = 1+ceil( g*(inimage[i]/lost_photons)  );
//              printf("  (%d,%d) auto-selected psf_radius = %d pixels\n",x0,y0,psf_radius);
            }
        }
        if(psf_radius == 0) psf_radius = 1;
        /* limit psf kernel to be no bigger than 4x the input image */
        maxwidth = xpixels;
        if(ypixels > maxwidth) maxwidth = ypixels;
        if(psf_radius > maxwidth) psf_radius = maxwidth;

        /* given the radius, how many photons will escape? */
        if(psftype == GAUSS)
        {
            /* r = sqrt(-log(lost_photons/total_photons)/log(16))*fwhm */
            /* lost_photons = total_photons*exp(-log(16)*(r^2/fwhm^2)) */
            rsq = psf_radius;
            rsq = rsq/fwhm_pixels;
            rsq = rsq*rsq;
            lost_photons = inimage[i]*exp(-log(16.)*rsq);
        }
        if(psftype == FIBER)
        {
            /* r ~ g*total_photons/lost_photons
               normalized integral from r=inf to "r" :  g/sqrt(g**2+r**2) */
            lost_photons = inimage[i]*g/sqrt(g*g+psf_radius*psf_radius);
        }
        /* accumulate this so we can add it to the whole image */
        total_lost_photons += lost_photons;

        for(dx=-psf_radius;dx<=psf_radius;++dx)
        {
            for(dy=-psf_radius;dy<=psf_radius;++dy)
            {
                /* this could be more efficient */
                k = kernel_size*(kernel_size/2+dy)+kernel_size/2+dx;
                if(kernel[k] == 0.0) continue;

                rsq = dx*dx+dy*dy;
                if(rsq > psf_radius*psf_radius) continue;
                x = x0+dx;
                y = y0+dy;
                if(x<0 || x>xpixels) continue;
                if(y<0 || y>ypixels) continue;

                /* index into output array */
                j = y*xpixels+x;
                /* do not wander off the output array */
                if(j<0 || j > pixels) continue;

                outimage[j] += inimage[i]*kernel[k];
            }
        }
    }
    /* now we have some lost photons, add them back "everywhere" */
    lost_photons = total_lost_photons/pixels;
    printf("adding back %g lost photons\n",total_lost_photons);
    for(i=0;i<pixels;++i)
    {
        outimage[i] += lost_photons;
    }

    /* don't need kernel anymore. but should we always allocate outimage? */
    free(kernel);
    return outimage;
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
            sq=sqrtf(2.0*xm);
            alxm=logf(xm);
            g=xm*alxm-gammln(xm+1.0);
        }
        do {
            do {
                /* y is a deviate from a lorentzian comparison function */
                y=tanf(M_PI*ran1(idum));
                /* shift and scale */
                em=sq*y+xm;
            } while (em < 0.0);         /* there are no negative Poisson deviates */
            /* round off to nearest integer */
            em=floor(em);
            /* ratio of Poisson distribution to comparison function */
            /* scale it back by 0.9 to make sure t is never > 1.0 */
            t=0.9*(1.0+y*y)*expf(em*alxm-gammln(em+1.0)-g);
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
        fac=sqrtf(-2.0*logf(rsq)/rsq);
        gset=v1*fac;
        iset=1;         /* we now have a spare deviate */
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
    return -logf(dum);
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

//          frame.mmapdata = mmap(NULL,2*frame.width*frame.height+frame.header_size,PROT_READ,MAP_SHARED,fileno(frame.handle),0);
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
