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

struct adjust *get_neighbours(double *npos, int lennd, int lentn, struct adjust *origin, double adrate);
void clear_ll(struct adjust *root);
void print_ll(struct adjust *root);

void som_train_loop(double *df, double *weights, double *distnd, Sint *prep, Sint *plendf,
		Sint *plennd, double *plrinit, double *freq, double *alpha, Sint *pdim, double *gt, double *npos, Sint *pradius,
		Sint *plentn, Sint *plentd, double *currtrain, Sint *plentr){

	int nearest, totiter, phase;
	int i, j, k, l, m, n, o, p;
	int w1, w2;
	double min, max;
	double meandist, adrate;
	struct adjust *nneigh, *nonneigh;
	double lrinit = *plrinit;
	int rep = *prep, lendf = *plendf, lennd = *plennd, dim = *pdim, initradius = *pradius;
	int lentn = *plentn, lentd = *plentd, lentr = *plentr;
	int nind;
	double dist, tmp, dm, lr, errorsum, radius;
	int nodegrow, x;
	struct adjust *root, *tp, *current, *tnode, *hptr, *hptr2, *newnode, *newnode_f;
	double sr = 0.5;

	root = (struct adjust *) malloc( sizeof(struct adjust) );
	root -> next = (struct adjust *) malloc( sizeof(struct adjust) );
	root -> next -> next = NULL;
	root -> next -> nodeid = -1;

	totiter = rep * lendf;

	phase = 1;

	// Loop over iterations of the input data
	for(i = 0; i<rep; i++){

		// Reset Frequencies
		// In order to be able to delete unneeded nodes
		for(j = 0; j<lennd; j++) freq[j] = 0;

		nodegrow = 0;
		errorsum = 0;

		// Reseting Learning rate during each iteration.
		// It should be considered if this:
		//  -Should be ommited for the training phase
		//  -The calculation of the learning rate from the traditional Kohonen should be used for
		//   the phase 2.
		lr = lrinit;

		// Loop over number of observations
		for(j = 0; j<lendf; j++){

			//Select Random observation
			x = (lendf*unif_rand());

			//Adjust learning Rate
			lr = *alpha * ( 1-(3.8/lennd))*lr;

			// Find best matching node
			nind = 0;
			nearest = -1;
			dm = INFINITY;

			for (k = 0; k < lennd; k++) {

				//Compute Euclidean Distance
				dist = 0.0;
				for (l = 0; l < dim; l++) {
					tmp = df[x + lendf*l] - weights[k + l*lentn];
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

			if(nearest == -1) error("Critical: No best matching node found. This should not happen.");

			//Update some counters
			distnd[nearest] += dm;
			errorsum += dm;
			freq[nearest]++;

			//Detect Radius
// This is broken because it gives back only the max radius and 1... Problem with integer and doubles.
			if(phase == 1) radius = initradius;
			else radius = initradius * ((rep - i) / rep) + 1.0;

			// Find neighbourhood
			root -> nodeid = nearest;
			root -> adrate = 1;
			root -> next = NULL;

			for(k = radius; k>0; k--){

				adrate = k/(radius+1.0);
				tnode = get_neighbours(npos, lennd, lentn, root, adrate);
				current = root;
				while(current -> next != NULL) current = current -> next;
				current -> next = tnode;

			}

			//Adjust neighbourhood
			current = root;
			while(current -> next != NULL){

				for(k; k<dim; k++){
					weights[current -> nodeid + k*lentn] = weights[current -> nodeid+ k*lentn] +
						(df[x + k*lendf] - weights[current -> nodeid+ k*lentn]) * lr * current -> adrate;
				}

				current = current -> next;
			}

			// Growth / Spreading
			if(distnd[nearest] > *gt && phase == 1){

				//Determine if maximum size of net is reached
				if(!(lentn-1 > lennd)) error("Number of nodes exceeded maximum capacity. Consider adjusting max gridsize or reducing Spreading Factor.");

				//Determine the number of direct neighbours.
				//This is kind of hacky...
				current = root;
				tmp = 0;
				while(current != NULL){
					if(current -> adrate >= (initradius/(initradius+1.0)) ) tmp++;
					current = current -> next;
				}

				current = root;

				if(tmp > 4){

					//Node has 4 direct neighbours. Growth is not possible.
					//Therefore the error is spread to neighbouring units.
					distnd[nearest] = distnd[nearest] / 2;
					for(l=1; l < 5; l++){

						current = current -> next;
						// Paper suggests values between 0 and 1
						distnd[current -> nodeid] = distnd[current -> nodeid] * (1 + sr);

					}

				} else {

					// Update counter for phase change.
					nodegrow += 1;

					//Check for growth on all 4 sides of nearest
					for(l=0; l<4; l++){

						//If tmp remains 0, an empty spot for a new node was found.
						tmp = 0;
						for(m=0; m < lennd; m++){
							if(npos[m] == npos[nearest] + ax[l] && npos[m + lentn] == npos[nearest + lentn] + ay[l]) tmp = 1;
						}

						//If true, generate new node on this position. If not skip to next possible position.
						if(tmp == 0){

							//Adjust number of nodes of GSOM
							lennd++;

							//Init values for new node.
							npos[lennd-1] = npos[nearest] + ax[l];
							npos[lennd-1 + lentn] = npos[nearest + lentn] + ay[l];
							distnd[lennd-1] = 0;
							freq[lennd - 1] = 0;

							//Get the neighbours of the new node.
							//struct adjust *newnode;
							newnode = malloc(sizeof(struct adjust));
							newnode -> next = NULL;
							newnode -> nodeid = lennd-1;

							nneigh = get_neighbours(npos, lennd, lentn, newnode, 0);

							clear_ll(newnode);

							//Set new weights for the new node. 4 Cases (A, B, C & D are considered)
							if(nneigh -> next != NULL){

								//Case B (More than one direct neighbour exists)
								for(n=0; n<dim; n++){
									weights[lennd-1 + n*lentn] = (weights[nneigh -> nodeid + n*lentn] + weights[nneigh -> next -> nodeid + n*lentn]) / 2;
								}

							} else{

								//Get the direct neighbours of nearest
								newnode_f = malloc(sizeof(struct adjust));
								newnode_f -> next = NULL;
								newnode_f -> nodeid = nearest;

								nonneigh = get_neighbours(npos, lennd, lentn, newnode_f, 0);
								clear_ll(newnode_f);

								//Sanity Check for debugging reasons
								if(nneigh-> nodeid != nearest) error("Topology of the GSOM is broken...");

								//Check the neighbours of nearest...
								hptr = nonneigh;
								tmp = -1;
								while(hptr != NULL){

									//...Remove new node from list
									if(hptr -> nodeid == lennd-1){ //nneigh -> nodeid
										hptr2 = hptr -> next;
										hptr -> nodeid = hptr2 -> nodeid;
										lr = lrinit;
										hptr -> next = hptr2 -> next;
										free(hptr2);
									}

									//...Determine if case A applies, tmp stores the relevant nodeid
									for(o=0; o < 4; o++){
										if(npos[hptr -> nodeid] == npos[lennd-1] + ax[o]*2 && npos[hptr -> nodeid + lentn] == npos[lennd-1 + lentn] + ay[o]*2 ) tmp = hptr -> nodeid;
									}

									hptr = hptr -> next;

								}

								w1 = nearest;
								if(tmp != -1){

									//Case A (Parent node w1 has a node w2 lying in the same direction as w1 lies in respect to new node)
									w2 = tmp;
									for(o=0; o < dim; o++){
										if(weights[w1 + o*lentn] < weights[w2 + o*lentn]){
											weights[lennd-1 + lentn*o] = weights[w1 + o*lentn]-(weights[w2 + o*lentn] - weights[w1 + o*lentn]);
										}else{
											weights[lennd-1 + lentn*o] = weights[w1 + o*lentn]-(weights[w2 + o*lentn] - weights[w1 + o*lentn]);
										}

									}

								}else if(nonneigh == NULL){

									//Case D (only direct neghbour of new nowde has no neighbours)
									//Initialize w. average weights according to paper.
									for(o=0; o<dim; o++){

										//Find min and max for each weight
										min= INFINITY;
										max= -INFINITY;
										for(p=0; p<lennd; p++){
											if(weights[p + o*lentn] > max) max = weights[p + o*lentn];
											if(weights[p + o*lentn] < min) min = weights[p + o*lentn];
										}

										weights[lennd-1 + lentn*o] = (min + max) / 2;

									}

								}else{

									//Case C (Parent node w1 has neighbours but none that qualifies for case A)
									w2 = nonneigh -> nodeid; //Random existing neighbour of nearest.

									for(o=0; o < dim; o++){
										if(weights[w1 + o*lentn] < weights[w2 + o*lentn]){
											weights[lennd-1 + lentn*o] = weights[w1 + o*lentn]-(weights[w2 + o*lentn] - weights[w1 + o*lentn]);
										}else{
											weights[lennd-1 + lentn*o] = weights[w1 + o*lentn]-(weights[w2 + o*lentn] - weights[w1 + o*lentn]);
										}

									}

								}

								clear_ll(nonneigh);

							}

							clear_ll(nneigh);

						}

					}

					//Reset Learningrate and Error of nearest
					lr = lrinit;
					distnd[nearest] = 0;
				}

				//End of Growth and distribution after this line.

			}

			//Element j has been checked. Reset LL.
			clear_ll(root);
			root = malloc(sizeof(struct adjust));

		}

		//Iteration over DF (well the no of observations anyways) is completed

		//Update Training Progress
		meandist = errorsum / lendf;
		currtrain[i + 0*lentr] = i;
		currtrain[i + 1*lentr] = phase;
		currtrain[i + 2*lentr] = meandist;
		currtrain[i + 3*lentr] = lennd;
		currtrain[i + 4*lentr] = nodegrow;

		//Remove Empty Units
		j=0;
		while(j < lennd && phase == 1){

			if(freq[j] > 0) j++;

			else{
				freq[j] = freq[lennd-1];
				freq[lennd-1] = 0;

				npos[j] = npos[lennd-1];
				npos[j + lentn] = npos[lennd-1 + lentn];
				npos[lennd-1] = 0.0/0.0;
				npos[lennd-1 + lentn] = 0.0/0.0;

				for(k=0; k<dim; k++){
					weights[j + k*lentn] = weights[lennd-1 + k*lentn];
					weights[lennd-1 + k*lentn] = 0.0/0.0;
				}

				lennd--;

			}

		}

		//Check if Phase should change
		if((i > 4 && currtrain[i-5 + 3*lentr] >= currtrain[i-1 + 3*lentr]) || (i+1) > (rep/2)){
			phase = 2;
		}

	}

	//Iteration i is completed

	//Update Return Values
	*plennd = lennd;

}


//This function searches for the neighbourhood of a node (required as struct adjust)
//Requires: Origin (struct adjust), Lenght of Nodes, Lenght of node Data Structure, Pointer to npos, and the adrate for the new points
//Duplicates (that already esixist in Origin) are sorted out and not returned. The returned LL only contains new elements, and a Pointer
struct adjust *get_neighbours(double *npos, int lennd, int lentn, struct adjust *origin, double adrate){

	struct adjust *nroot, *tmp;
	int isneighbour, exclude;
	int l, m;

	nroot = NULL;

	//Origin is a linked list containing of n nodes, to which neighbours should be found
	if(origin == NULL) error("Linked list which should contain the origin is empty.");

	while(origin != NULL){

		//Check all nodes are to be added to the new ll
		for(l=0; l<lennd; l++){

			isneighbour = 0;

			//Check if the element is a neighbour
			for(m=0; m < 4; m++){

				if(npos[l] == npos[origin -> nodeid] + ax[m] && npos[l + lentn] == npos[origin -> nodeid + lentn] + ay[m])
					isneighbour = 1;

			}

			if(isneighbour == 1){

				//Sort out nodes that are in the input LL
				exclude=0;
				tmp = origin;
				while(tmp -> next != NULL){
					if(tmp -> nodeid == l) exclude = 1;
					tmp = tmp -> next;
				}

				//Add Node to new LL
				if(exclude == 0){
					//printf("Added node %d to linked list.\n", l);
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

	return nroot;

}


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

//For debugging
void print_ll(struct adjust *root){
	while(root != NULL){
		printf("NodeId = %d, Adrate=%f", root -> nodeid, root -> adrate);
		root = root -> next;
	}
}
