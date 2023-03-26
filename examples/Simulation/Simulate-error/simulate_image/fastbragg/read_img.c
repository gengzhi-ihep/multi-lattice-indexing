#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

   fp = fopen("new.txt", "w");
   int n;
   for(n=0; n<frame_width*frame_height; n++){
      fprintf(fp, "%d\n", data[n]);
   }
   fclose(fp);

   free(header);
   free(data);

   return;

}
