////////////////////////////////////////
//trainloop.c - part of GSOM.r
//Alex Hunziker - 17.10.2016
////////////////////////////////////////

#include <R.h>

#define EPS 1e-4                /* relative test of equality of distances */

struct adjust{
	int nodeid;
	double adrate;
	struct adjust *next;
};

void som_train_loop(double *df, double *weights, double *distnd, Sint *prep, Sint *plendf
		Sint *plennd, Sint *plrinit, double *freq, double *alhpa, Sind *pdim, Sind *pgt){
	
	int winner, totiter;
	int nodegrow, errorsum, x, tmp, dm;
	int to, bo, le, ri;
	int rep = prep, lendf = plendf, lrinit = plr, lennd = plennd, dim = pdim, gt = pgt;
	struct adjust *root, *tp, *current;
	
	root = (struct adjust *) malloc( sizeof(struct adjust) ); 
	root -> next = NULL;
	
	totiter = rep * lendf;
	
	// Loop over iterations of the input data
	for(i = 1, i<=rep, i++){
		
		// Reset Frequencies 
		// In order to be able to delete unneeded nodes
		for(j = 0, j<lennd, j++) freq[j] = 0;
		
		nodegrow = 0;
		errorsum = 0;
		
		lr = lrinit;
		
		// Loop over number of observations
		for(j = 0, j<lendf, j++){
			
			x = (int)(lendf * unif_rand());
			
			lr = alpha * ( 1-( 3.8/lennd ) ) * lr;
			
			// Find best matching node
			nind = 0;
			winner = -1;
			dm = DOUBLE_XMAX
			for (k = 0; k < lennd; k++) {
				dist = 0.0;
				for (l = 0; l < dim; l++) {
					tmp = data[x + lendf*l] - weights[x + l*lennd];
					dist += tmp * tmp;
				}

				if (dist <= dm * (1 + EPS)) {
					if (dist < dm * (1 - EPS)) {
						nind = 0;
						nearest = cd;
					} else {
						if(++nind * unif_rand() < 1.0) nearest = cd;
					}
					dm = dist;
				}
			}
			
			distnd[nearest] += dm;
			errorsum += errorsum;
			freq[nearest]++;
			
			radius = "MISSING FEATURE";
			
			// Find neighbourhood
			for(k; k>0, k--){
				current = root;
				while(current->next != NULL){
					for(l=1; l<lennd; l++){
						if(npos[l] == npos[current -> nodeid]+1 && npos[l+1] == npos[current -> nodeid+1] ||
							npos[l] == npos[current -> nodeid]-1 && npos[l+1] == npos[current -> nodeid+1] ||
							npos[l] == npos[current -> nodeid] && npos[l+1] == npos[current -> nodeid+1]+1 ||
							npos[l] == npos[current -> nodeid] && npos[l+1] == npos[current -> nodeid+1]-1){
							
							temp=0;
							tp = root;
							while(tp != NULL){
								if(tp -> nodeid == l) temp = 1;
								tp = tp -> next;
							}
							if(temp == 0){
								current -> next = (struct adjust *) malloc( sizeof(struct adjust) );
								current -> nodeid = l;
								current -> adrate = 1 - l/5;
								current = current -> next;
							}
						}
					}
				}
			}
			
			//Adjust neighbourhood
			while(root -> next != NULL){
				for(k; k<dim; k++){
					weights[nearest+ k*lennd] = weights[nearest+ k*lennd] + 
						(data[x + k*lendf] - weights[nearest+ k*lennd]) * lr;
				}
				tp=root;
				root = root -> next;
				free(tp);
			}
			
			// Growth / Spreading
			if(distnd[nearest] > gt){
				current = root;
				tmp = 0;
				while(current -> adrate >= 1 - 2/5){
					tmp ++;
				}
				if(tmp > 4){
					distnd[nearest] = distnd[nearest] / 2;
					for(l=1; l < 5; l++){
						// Paper suggests values between 0 and 1
						current = current -> next;
						distnd[current -> nodeid] = distnd[current -> nodeid] * (1 + 0.5);
					}
				} else {
					nodegrow += 1;
					
					//Growthcondition
					
					lr = lrinit;
					dist[nearest] = 0;
					
				}
			}
			
		}
		
	}
}