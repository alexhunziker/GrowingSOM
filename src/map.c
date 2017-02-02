////////////////////////////////////////
//map.c - GrowingSOM
//Alex Hunziker - 2017
////////////////////////////////////////

#include <R.h>
#include <math.h>

#define EPS 1e-4                /* relative test of equality of distances */

void map_data(Sint *plendf, Sint *lennd, Sint *dim, double *df, double *codes, double *bmn, double *ndist, double *freq){

  int lendf = *plendf;

  //Declare Variables
  int i, k, l, nind, nearest;
  double dm, dist, tmp;

  //Loop through observations to match
  for(i=0; i<lendf; i++){

    nind = 0;
    nearest = -1;
    dm = INFINITY;

    //Check every node
    for (k = 0; k < *lennd; k++) {

      //Compute Euclidean Distance
      dist = 0.0;
      for (l = 0; l < *dim; l++) {
        tmp = df[i + lendf*l] - codes[k + l * *lennd];
        dist += tmp * tmp;
      }

      //Check if current node is nearest so far
      if (dist <= dm * (1 + EPS)) {
        if (dist < dm * (1 - EPS)) {
          nind = 0;
          nearest = k;
        } else {
          if(++nind * unif_rand() < 1.0) nearest = k;
        }
        dm = dist;
      }

    }

    ndist[i] = sqrt(dm / (double)*dim);
    bmn[i] = nearest+1;

    if(nearest == -1) error("No nearest unit found.");

    freq[nearest] = freq[nearest] + 1;

  }
}
