
#include <stdlib.h> /* malloc, free */
#include <stdio.h> /* printf, fprintf, feof, fgets */
#include <libgen.h> /* basename */
#include <string.h> /* strlen */
#include <ctype.h> /* isspace */

#define MAXCHAR 1024
int main(int argc, char** argv) {
   int ret = 0;
   int cnt;
   char data[MAXCHAR];
   double* nodes = NULL; int n_nodes;
   int* elems = NULL; int n_elems;

   if(argc != 2) { fprintf(stderr,"usage: %s <netgen.vol>\n",basename(argv[0])); return 1; }
   FILE* F = fopen(argv[1],"r");
   if(!F) { fprintf(stderr,"error: failed to open %s\n",argv[0]); goto __quit; }
   printf("reading %s\n",argv[1]);

   if(!fgets(data, 8, F)) goto __quit;
   /* DANGER: inconsistent use of \r and \n by netgen */
   if(!strcmp(data,"mesh3d\r")) {fprintf(stderr,"err: bad header\n"); goto __quit;}
   while(!feof(F)) {
      if(!fgets(data, MAXCHAR, F)) break;

      if(!strcmp(data,"dimension\n")) {
         int dim = 0;
         cnt = fscanf(F,"%d\n",&dim);
         if((cnt != 1) || (dim != 3)) {fprintf(stderr,"err: bad dimension\n"); goto __quit;}
         printf("dimension %d\n",dim);
      }
      else if(!strcmp(data,"geomtype\n")) {
         int type = -1;
         cnt = fscanf(F,"%d\n",&type);
         if((cnt != 1) || (type != 0)) {fprintf(stderr,"err: bad geomtype\n"); goto __quit;}
         printf("geomtype %d\n",type);
      }
      else if(!strcmp(data,"points\n")) {
         cnt = fscanf(F,"%d\n",&n_nodes);
         printf("%d points\n",n_nodes);
         int node = 0;
         nodes = malloc(sizeof(double)*n_nodes*3);
         if(!nodes) return -1;
         while(!feof(F) && node < n_nodes ) {
            cnt += fscanf(F," %lf %lf %lf\n", &nodes[3*node+0], &nodes[3*node+1], &nodes[3*node+2]);
            node++;
         }
         if(cnt != n_nodes*3+1) {fprintf(stderr,"err: bad points\n"); goto __quit;}
      }
      else if(!strcmp(data,"volumeelements\n")) {
         cnt = fscanf(F,"%d\n",&n_elems);
         printf("%d volumeelements\n",n_elems);
         int elem = 0;
         elems = malloc(sizeof(double)*n_elems*4);
         if(!elems) return -1;
         while(!feof(F) && elem < n_elems ) {
            cnt += fscanf(F,"1 4 %d %d %d %d\n", &elems[4*elem+0], &elems[4*elem+1], &elems[4*elem+2], &elems[4*elem+3]);
            elem++;
         }
         if(cnt != n_elems*4+1) {fprintf(stderr,"err: bad elems\n"); goto __quit;}
      }
   }
//   if((n_elems > 0) && (n_nodes > 0)) {
//      printf("happiness!\n");
//   }

__quit:
   free(nodes); nodes = NULL;
   free(elems); elems = NULL;
   if(F) {
      ret = fclose(F);
      if(ret) { fprintf(stderr,"error: failed to close %s\n",argv[1]); }
   }
}
