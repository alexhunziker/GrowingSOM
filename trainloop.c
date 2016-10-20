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

*struct adjust get_neighbours(double *npos, struct adjust *origin);
void clear_ll(struct adjust *root);

void som_train_loop(double *df, double *weights, double *distnd, Sint *prep, Sint *plendf,
		Sint *plennd, double *plrinit, double *freq, double *alpha, Sint *pdim, double *gt, double *npos, Sint *pradius){

	//General Variables
	int nearest, totiter, phase;

	//Variables for Iterations
	int i, j, k, l, m;

	//Variables For Nodegrowth


	//Variables From Pointers
	int rep = *prep, lendf = *plendf, lrinit = *plrinit, lennd = *plennd, dim = *pdim, initradius = *pradius;

	//Variables Fror Neighbour Detection
	int ax[4] = {0, 0, 1, -1};
	int ay[4] = {1, -1, 0, 0};

	int nind, radius;
	double dist, tmp, dm, lr, errorsum;
	int nodegrow, x;

	struct adjust *root, *tp, *current, *tnode, *hptr, *hptr2;

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
				while(current -> next != NULL){
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
						tmp = 0;
						for(m=0; m < lennd; m++){
							if(npos[m] == npos[nearest + ax[l]] && npos[m + lennd] == npos[nearest + ay[l] + lennd]) tmp = 1;
						}
						if(tmp == 0){
							lennd++;
							npos[lennd-1] = nearest + ax[l];
							npos[lennd*2 - 1] = nearest + ay[l];
							distnd[lennd-1] = 0;
							freq[lennd - 1] = 0;

							//Get the neighbours of the new node.
							struct adjust *newnode;
							newnode -> next = NULL;
							newnode -> nodeid = m;
							nneigh = get_neighbours(npos, lennd, newnode);
							clear_ll(newnode);

							if(nneight -> next != NULL){
								//Case B
								for(n=0; n<dim; n++){
									nweight[n] = weights[npos[nneigh] + lennd*n] + weights[npos[nneigh + lennd] + lennd*n] / 2;
								}

							} else{
								//Get the neighbours of the First Neighbour of New node
								struct adjust *newnode_f;
								newnode_f -> next = NULL;
								newnode_f -> nodeid = nneigh -> nodeid;
								nonneigh = get_neighbours(npos, lennd, newnode_f);
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
										if(weights[w1 + o*lennd] < weights[w2 + o*lennd]){
											weights[lennd-1 + lennd*o] = weights[w1 + o*lennd]-(weights[w2 + o*lennd] - weights[w1 + o*lennd]);
										}else{
											weights[lennd-1 + lennd*o] = weights[w1 + o*lennd]-(weights[w2 + o*lennd] - weights[w1 + o*lennd]);
										}
									}

								}else if(nonneigh == NULL){
									//Case D
									for(o=0; o<dim; o++){
										min=infinity();
										max=-infinity();
										for(p=0; p<lennd; p++){
											if(weights[p + o*lennd] > max) max = weights[p + o*lennd];
											if(weights[p + o*lennd] < min) min = weights[p + o*lennd];
										}
										weights[lennd-1 + lennd*o] = (min + max) / 2;
									}

								}else{

									//Case C
									w2 = nonneigh -> nodeid;
									for(o=0, o < dim; o++){
										if(weights[w1 + o*lennd] < weights[w2 + o*lennd]){
											weights[lennd-1 + lennd*o] = weights[w1 + o*lennd]-(weights[w2 + o*lennd] - weights[w1 + o*lennd]);
										}else{
											weights[lennd-1 + lennd*o] = weights[w1 + o*lennd]-(weights[w2 + o*lennd] - weights[w1 + o*lennd]);
										}
									}
								}
								clear_ll(nonneigh);
							}
							clear_ll(nneigh);
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
		while(j < lennd){
			if(freq[j] > 0) j++;
			else{
				freq[j] = freq[lennd-1];
				freq[lennd-1] = na();

				npos[j] = npos[lennd-1];
				npos[j + lennd] = npos[lennd-1 + lennd];
				npos[lennd-1] = na();
				npos[lennda-1 + lennd] = na();

				for(k=0; k<dim; k++){
					weights[j + k*lennd] = weights[lennd + k*lennd];
					weights[lennd + k*lennd]; = na();
				}

				lennd--;
			}
		}

		//Check if Phase should change

	}
}

struct *adjust get_neighbours(double *npos, int lennd, struct adjust *origin, int adrate=0){

	struct adjust *nroot, *tmp;
	int isneighbour, exclude;
	int l, m;

	nroot = NULL;

	//Origin is a list of n nodes, to which neighbours should be found
	while(orign != NULL){
		for(l=0; l<lennd; l++){
			isneighbour = 0;

			for(m=0; m < 4; m++){
				if(npos[l] == npos[nearest + ax[m]] && npos[l + lennd] == npos[nearest + ay[m] + lennd])
					isneighbour = 1;
			}

			if(isneighbour == 1){
				printf("Found a neighbour\n");

				//Sort out nodes that are in the input LL
				exclude=0;
				tmp = root;
				while(tmp -> next != NULL){
					if(tmp -> nodeid == l) exclude = 1;
					tmp = tmp -> next;
				}

				//Add Node to new LL
				if(exclude == 0){
					tmp = (struct adjust *) malloc( sizeof(struct adjust) );
					tmp -> next = nroot;
					tnode -> nodeid = l;
					tnode -> adrate = adrate;
					nroot = tmp;
				}

			}
		}
		origin = origin -> next;
	}

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

void clear_ll(struct adjust *root){
	struct adjust *tmp;
	while(root != NULL){
		tmp = root;
		root = root -> next;
		free(tmp);
	}
	return;
}
