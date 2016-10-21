////////////////////////////////////////
//trainloop.c - part of GSOM.r
//Alex Hunziker - 17.10.2016
////////////////////////////////////////

#include <R.h>
#include <unistd.h>
#include <math.h>

#define EPS 1e-4                /* relative test of equality of distances */

struct adjust{
	int nodeid;
	double adrate;
	struct adjust *next;
};

int ax[4] = {0, 0, 1, -1};
int ay[4] = {1, -1, 0, 0};

struct adjust *get_neighbours(double *npos, int lennd, int lentn, struct adjust *origin, int adrate);
void clear_ll(struct adjust *root);

void som_train_loop(double *df, double *weights, double *distnd, Sint *prep, Sint *plendf,
		Sint *plennd, double *plrinit, double *freq, double *alpha, Sint *pdim, double *gt, double *npos, Sint *pradius,
		Sint *plentn, Sint *plentd, double *currtrain, Sint *plentr){

	int nearest, totiter, phase;
	int i, j, k, l, m, n, o, p;
	int w1, w2;
	double min, max;
	double meandist;
	struct adjust *nneigh, *nonneigh;
	int rep = *prep, lendf = *plendf, lrinit = *plrinit, lennd = *plennd, dim = *pdim, initradius = *pradius;
	int lentn = *plentn, lentd = *plentd, lentr = *plentr;
	int nind, radius;
	double dist, tmp, dm, lr, errorsum;
	int nodegrow, x;
	struct adjust *root, *tp, *current, *tnode, *hptr, *hptr2, *newnode, *newnode_f;

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
			dm = INFINITY;
			printf("lennd %d", lennd);
			for (k = 0; k < lennd; k++) {
				printf("l");
				dist = 0.0;
				for (l = 0; l < dim; l++) {
					tmp = df[x + lendf*l] - weights[k + l*lentn];
					dist += tmp * tmp;
				}

				if (dist <= dm * (1 + EPS)) {
					printf("i");
					if (dist < dm * (1 - EPS)) {
						printf("il");
						nind = 0;
						nearest = k;
					} else {
						if(++nind * unif_rand() < 1.0) nearest = k;
					}
					dm = dist;
				}
			}
			if(nearest == -1) error("Critical: No best matching unit found.");

			distnd[nearest] += dm;
			errorsum += dm;
			freq[nearest]++;

			//Detect Radius.
			if(phase == 1) radius = initradius;
			else radius = initradius * ((rep - i) / rep) + 1;

			// Find neighbourhood
			root -> nodeid = nearest;
			root -> adrate = 1;
			root -> next = NULL;

			for(k = radius; k>0; k--){

				tnode = get_neighbours(npos, lennd, lentn, root, 1); // Fix adrate.
				current = root;
				while(current -> next != NULL) current = current -> next;
				current -> next = tnode;

			}

			//Adjust neighbourhood
			current = root;
			while(current -> next -> next != NULL){
				for(k; k<dim; k++){
					weights[current -> nodeid + k*lentn] = weights[current -> nodeid+ k*lentn] +
						(df[x + k*lendf] - weights[current -> nodeid+ k*lentn]) * lr * current -> adrate;
				}
				current = current -> next;
			}

			// Growth / Spreading
			if(distnd[nearest] > *gt && phase == 1){
				current = root;
				tmp = 0;
				// This is Broken.
				while(root -> next -> next != NULL){
					if(root -> adrate >= 1 - 2/5) tmp ++;
					tp=root;
					root = root -> next;
					free(tp);
				}

				if(tmp == 4){
					//Node has 4 direct neighbours. Growth is not possible.
					//Therefore the error is spread to neighbouring units.
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
						tmp = 0;
						for(m=0; m < lennd; m++){
							if(npos[m] == npos[nearest] + ax[l] && npos[m + lentn] == npos[nearest + lentn] + ay[l]) tmp = 1;
							printf("npos[m]: %f, npos[m + lentn]: %f\n", npos[m], npos[m + lentn]);
						}
						printf("tmp is: %f", tmp);
						if(tmp == 0){

							//Init New Node Values except weights
							lennd++;
							npos[lennd-1] = npos[nearest] + ax[l];
							printf("NEWNEWNEW %d \n", nearest + ax[l]);
							npos[lennd-1 + lentn] = npos[nearest + lentn] + ay[l];
							distnd[lennd-1] = 0;
							freq[lennd - 1] = 0;

							//Get the neighbours of the new node.
							//struct adjust *newnode;
							newnode = malloc(sizeof(struct adjust));
							newnode -> next = NULL;
							newnode -> nodeid = m;
							printf("newnode: %p, %d", newnode, newnode -> nodeid);
							printf("START INTERESSTING PART \n");
							nneigh = get_neighbours(npos, lennd, lentn, newnode, 0);
							clear_ll(newnode);

							//printf("Prt: %d", nneigh -> nodeid);
							if(NULL==NULL) printf("nneigh: is not null");
							else printf("It's null. That shouldn't happen.");

							if(nneigh -> next != NULL){
								//Case B
								for(n=0; n<dim; n++){
									weights[lennd-1 + n*dim] = (weights[nneigh -> nodeid + n*dim] + weights[nneigh -> next -> nodeid + n*dim]) / 2;
								}

							} else{
								//Get the neighbours of the First Neighbour of New node
								newnode_f = malloc(sizeof(struct adjust));
								newnode_f -> next = NULL;
								newnode_f -> nodeid = nneigh -> nodeid;
								nonneigh = get_neighbours(npos, lennd, lentn, newnode_f, 0);
								clear_ll(newnode_f);

								//Delete New Node from Neighbours of New Nodes neighbours
								// & determine if case A applies.
								hptr = nonneigh;
								tmp = -1;
								while(hptr != NULL){
									if(hptr -> nodeid == nneigh -> nodeid){
										hptr2 = hptr -> next;
										hptr -> nodeid = hptr2 -> nodeid;
										hptr -> next = hptr2 -> next;
										free(hptr2);
									}
									for(o=0; o < 4; o++){
										if(npos[hptr -> nodeid] == npos[nearest] + ax[o]*2 && npos[hptr -> nodeid + lennd] == npos[nearest + lennd] + ay[o]*2 ) tmp = o;
									}
									hptr = hptr -> next;
								}

								w1 = nneigh -> nodeid;
								if(tmp != -1){
									//Case A
									hptr = nonneigh;
									for(o=0; o<tmp; o++){
										hptr = hptr -> next;
									}
									w2 = hptr -> nodeid;
									for(o=0; o < dim; o++){
										if(weights[w1 + o*lentn] < weights[w2 + o*lentn]){
											weights[lennd-1 + lennd*o] = weights[w1 + o*lentn]-(weights[w2 + o*lentn] - weights[w1 + o*lentn]);
										}else{
											weights[lennd-1 + lennd*o] = weights[w1 + o*lentn]-(weights[w2 + o*lentn] - weights[w1 + o*lentn]);
										}
									}

								}else if(nonneigh == NULL){
									//Case D
									for(o=0; o<dim; o++){
										min= INFINITY;
										max= -INFINITY;
										for(p=0; p<lennd; p++){
											if(weights[p + o*lentn] > max) max = weights[p + o*lentn];
											if(weights[p + o*lentn] < min) min = weights[p + o*lentn];
										}
										weights[lennd-1 + lennd*o] = (min + max) / 2;
									}

								}else{

									//Case C
									w2 = nonneigh -> nodeid;
									for(o=0; o < dim; o++){
										if(weights[w1 + o*lentn] < weights[w2 + o*lentn]){
											weights[lennd-1 + lennd*o] = weights[w1 + o*lentn]-(weights[w2 + o*lentn] - weights[w1 + o*lentn]);
										}else{
											weights[lennd-1 + lennd*o] = weights[w1 + o*lentn]-(weights[w2 + o*lentn] - weights[w1 + o*lentn]);
										}
									}
								}
								clear_ll(nonneigh);
							}
							clear_ll(nneigh);
							error("Yep we are growing a new one.");
						}
					}

					lr = lrinit;
					distnd[nearest] = 0;
				}
			}

		}

		//Update Training Progress
		meandist = errorsum / lendf;
		currtrain[i + 0*lentr] = i;
		currtrain[i + 1*lentr] = phase;
		currtrain[i + 2*lentr] = meandist;
		currtrain[i + 3*lentr] = lennd;
		currtrain[i + 4*lentr] = nodegrow;

		//Remove Empty Units
		while(j < lennd && phase == 1){
			if(freq[j] > 0) j++;
			else{
				freq[j] = freq[lennd-1];
				freq[lennd-1] = 0;

				npos[j] = npos[lennd-1];
				npos[j + lennd] = npos[lennd-1 + lennd];
				npos[lennd-1] = 0.0/0.0;
				npos[lennd-1 + lennd] = 0.0/0.0;

				for(k=0; k<dim; k++){
					weights[j + k*lentn] = weights[lennd + k*lentn];
					weights[lennd + k*lentn] = 0.0/0.0;
				}

				lennd--;
			}
		}

		//Check if Phase should change
		if(i > 4 && currtrain[i-5 + 3*i] >= currtrain[i + 3*i] || i > (rep/2)){
			phase = 2;
		}

	}
}

struct adjust *get_neighbours(double *npos, int lennd, int lentn, struct adjust *origin, int adrate){

	struct adjust *nroot, *tmp;
	int isneighbour, exclude;
	int l, m;

	printf("Search for n of: %d", origin -> nodeid);

	nroot = NULL;

	//Origin is a list of n nodes, to which neighbours should be found
	if(origin == NULL) error("Origin is null.");
	while(origin != NULL){
		for(l=0; l<lennd; l++){
			printf("loop::%f::", npos[origin -> nodeid]);
			isneighbour = 0;

			for(m=0; m < 4; m++){
				if(npos[l] == npos[origin -> nodeid] + ax[m] && npos[l + lentn] == npos[origin -> nodeid + lentn] + ay[m])
					{isneighbour = 1; printf("Candidate Found.");}
					printf("_cc_");
					printf("%f, %f", npos[origin -> nodeid] + ax[m], npos[origin -> nodeid + lentn] + ay[m]);
			}


			if(isneighbour == 1){
printf("_aa_");
				//Sort out nodes that are in the input LL
				exclude=0;
				tmp = origin;
				while(tmp -> next != NULL){
					if(tmp -> nodeid == l) exclude = 1;
					printf("Exclude Check, for found element.");
					tmp = tmp -> next;
				}
printf("_dd_%d", exclude);
				//Add Node to new LL
				if(exclude == 0){
					printf("Added to LL");
					tmp = (struct adjust *) malloc( sizeof(struct adjust) );
					tmp -> next = nroot;
					tmp -> nodeid = l;
					tmp -> adrate = adrate;
					nroot = tmp;
				}

			}
		}
		origin = origin -> next;

	}
	printf("nroot: %p", nroot);
	return nroot;
}

	/* Old Neighbour Function
	//For all Elements Provided
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
} */

// House keeping, to avoid memory leaks.
void clear_ll(struct adjust *root){
	struct adjust *tmp;
	while(root != NULL){
		tmp = root;
		root = root -> next;
		free(tmp);
	}
	return;
}
