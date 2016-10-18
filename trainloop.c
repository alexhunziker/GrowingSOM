////////////////////////////////////////
//trainloop.c - part of GSOM.r
//Alex Hunziker - 17.10.2016
////////////////////////////////////////

#include <R.h>
#include <unistd.h>

#define EPS 1e-4                /* relative test of equality of distances */

struct adjust{
	int nodeid;
	double adrate;
	struct adjust *next;
};

void som_train_loop(double *df, double *weights, double *distnd, Sint *prep, Sint *plendf,
		Sint *plennd, double *plrinit, double *freq, double *alpha, Sint *pdim, double *gt, double *npos, Sint *pradius){
	
	int nearest, totiter, nind, phase, radius;
	double dist, tmp, dm, lr, errorsum;
	int nodegrow, x;
	int i, j, k, l, m;
	int ax[4] = {0, 0, 1, -1};
	int ay[4] = {1, -1, 0, 0};
	int rep = *prep, lendf = *plendf, lrinit = *plrinit, lennd = *plennd, dim = *pdim, initradius = *pradius;
	struct adjust *root, *tp, *current, *tnode;
	
	root = (struct adjust *) malloc( sizeof(struct adjust) ); 
	root -> next = (struct adjust *) malloc( sizeof(struct adjust) );
	root -> next -> next = NULL;
	root -> next -> nodeid = -1; 
	
	totiter = rep * lendf;

	phase = 1;
	
	// Loop over iterations of the input data
	for(i = 1; i<=rep; i++){
		
		// Reset Frequencies 
		// In order to be able to delete unneeded nodes
		for(j = 0; j<lennd; j++) freq[j] = 0;
		
		nodegrow = 0;
		errorsum = 0;
		
		lr = lrinit;
		
		// Loop over number of observations
		for(j = 0; j<lendf; j++){
			
			//x = (lendf*unif_rand());
			x=4;
			
			lr = *alpha * ( 1-(3.8/lennd))*lr;
			
			// Find best matching node
			nind = 0;
			nearest = -1;
			dm = 9999;
			for (k = 0; k < lennd; k++) {
				dist = 0.0;
				for (l = 0; l < dim; l++) {
					tmp = df[x + lendf*l] - weights[k + l*lennd];
					dist += tmp * tmp;
				}

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
			
			distnd[nearest] += dm;
			errorsum += dm;
			freq[nearest]++;

			//Detect Radius.
			if(phase == 1) radius = initradius;
			else radius = initradius * ((rep - i) / rep) + 1;
			
			// Find neighbourhood
			root -> nodeid = nearest;
			root -> adrate = 1;
			for(k = radius; k>0; k--){
				current = root;
				if(root -> next == NULL) error("Error in Linked List");
				
				printf("Current: %d, %d\n", current -> nodeid, k);
				while(current->next != NULL){
					for(l=0; l<lennd; l++){
						printf("nops(l) %f, (l+1) %f", npos[l], npos[l+lennd]);
						if(npos[l] == npos[current -> nodeid]+1 && npos[l+lennd] == npos[current -> nodeid+lennd] ||
							npos[l] == npos[current -> nodeid]-1 && npos[l+lennd] == npos[current -> nodeid+lennd] ||
							npos[l] == npos[current -> nodeid] && npos[l+lennd] == npos[current -> nodeid+lennd]+1 ||
							npos[l] == npos[current -> nodeid] && npos[l+lennd] == npos[current -> nodeid+lennd]-1){
							
							printf("Found a neighbour\n");							
							tmp=0;
							tp = root;
							while(tp -> next != NULL){
								if(tp -> nodeid == l) tmp = 1;
								tp = tp -> next;
							}

							if(tmp == 0){
								tnode = (struct adjust *) malloc( sizeof(struct adjust) );
								tnode -> next = root;
								tnode -> nodeid = l;
								tnode -> adrate = 1 - l/5;
								root = tnode;
							}
						}
					}
					current = current -> next;
					if(root -> next == NULL) error("Error in Linked List");
				}
				printf("Is.null\n");
			}
			
			//Adjust neighbourhood
			current = root;
			while(current -> next -> next != NULL){
				for(k; k<dim; k++){
					weights[current -> nodeid + k*lennd] = weights[current -> nodeid+ k*lennd] + 
						(df[x + k*lendf] - weights[current -> nodeid+ k*lennd]) * lr * current -> adrate;
				}
				current = current -> next;
			}
			
			// Growth / Spreading
			if(distnd[nearest] > *gt && phase == 1){
				current = root;
				tmp = 0;
				while(root -> next -> next != NULL){
					if(root -> adrate >= 1 - 2/5) tmp ++;
					tp=root;
					root = root -> next;
					free(tp);
				}

				if(tmp == 4){
					printf("Dist\n");
					distnd[nearest] = distnd[nearest] / 2;
					for(l=1; l < 5; l++){
						// Paper suggests values between 0 and 1
						current = current -> next;
						distnd[current -> nodeid] = distnd[current -> nodeid] * (1 + 0.5);
					}
				} else {
					printf("Grow\n");
					nodegrow += 1;
					
					//Growthcondition
					for(l=0; l<4; l++){
						/*tmp = 0;
						for(m=0; m < lennd; m++){
							if(npos[m] == npos[nearest + ax[l]] && npos[m + lennd] == npos[nearest + ay[l] + lennd]) tmp = 1;
						}
						if(tmp == 0){
							lennd++;
							npos[lennd-1] = nearest + ax[l];
							npos[lennd*2 - 1] = nearest + ay[l];
						}*/
					}
					
					lr = lrinit;
					distnd[nearest] = 0;
					
				}
			}
			
		}
		
	}
}
