#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

double ValueOf(const char *keyword, char *header){
    double value;
    char *string;
    int keylen = strlen(keyword);

    string = header;

    while((char*) strstr(string, keyword) != NULL){
        string = (char *) strstr(string,keyword) + keylen;
    }

    string = (char *) strstr(string, "=");
    if(string == NULL) return 0.0;
    ++string;

    value = atof(string);

    return value;

}

void main(){

   char *filename = "noiseimage.img";

   char *header = (char *) malloc(1024*sizeof(char));

   FILE * fp = fopen(filename, "rb");

   fread(header, 512, 1, fp);

   int header_size = (int) ValueOf("HEADER_BYTES", header);
   int frame_width = (int) ValueOf("SIZE1", header);
   int frame_height = (int) ValueOf("SIZE2", header);

   unsigned short int *data = (unsigned short int *) malloc(sizeof(unsigned short int)*frame_width*frame_height);

   fseek(fp, 0, SEEK_CUR);

   fread(data, sizeof(unsigned short int), frame_width*frame_height, fp);

   fclose(fp);

   fp = fopen("new.txt", "r");
   int n;
   for(n=0; n<frame_width*frame_height; n++){
      fscanf(fp, "%d\n", &data[n]);
   }
   fclose(fp);

   hid_t file_id, dset_id, dspace_id;
   hsize_t dims[2] = {frame_width, frame_height};
   herr_t status;

   file_id = H5Fcreate("new.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   dspace_id = H5Screate_simple(2, dims, NULL);
   dset_id = H5Dcreate2(file_id, "data", H5T_NATIVE_USHORT, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   status = H5Dwrite(dset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

   H5Dclose(dset_id);
   H5Sclose(dspace_id);
   H5Fclose(file_id);

   free(header);
   free(data);

   return;

}
