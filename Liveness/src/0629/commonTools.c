/***************************************************************************
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "commonTools.h"

/* init feature vectors */
void initFeatures(TypeFeatureList *features) {
	features->nbFrame = features->vectDim = 0;
	features->featureVectors = NULL;
}

/* create feature vectors and allocate memory */
int createFeatures(TypeFeatureList *features, long vectDim, long nbFrame) {
	long i;

	features->nbFrame = nbFrame;
	features->vectDim = vectDim;
	features->featureVectors = (double **) malloc (sizeof(double *) * features->nbFrame);
	if (features->featureVectors == NULL) {
		fprintf(stderr, "ERROR: Couldn't allocate memory for data buffer (features->featureVectors)\n");
		exit(EXIT_FAILURE);
	}
	for(i=0; i<features->nbFrame; i++) {
		features->featureVectors[i] = (double *) malloc (sizeof(double) * features->vectDim);
		if (features->featureVectors[i] == NULL) {
			fprintf(stderr, "ERROR: Couldn't allocate memory for data buffer (features->featureVectors[%ld])\n", i);
			exit(EXIT_FAILURE);
		}
	}
	return 0;
}

/* reset feature vectors and free allocated memory */
void deleteFeatures(TypeFeatureList *features) {
	features->nbFrame = features->vectDim = 0;
	if (features->featureVectors != NULL) {
		free(features->featureVectors);
		features->featureVectors = NULL;
	}
}

/* copy feature vectors */
int copyFeatures(TypeFeatureList *source_feats, TypeFeatureList *target_feats) {
	long i, d;

	if ((target_feats->featureVectors != NULL) && (target_feats->nbFrame != 0))
		deleteFeatures(target_feats);
	createFeatures(target_feats, source_feats->vectDim, source_feats->nbFrame);

	for(i=0; i<source_feats->nbFrame; i++) {
	        for (d=0; d<source_feats->vectDim; d++) 
			target_feats->featureVectors[i][d] = source_feats->featureVectors[i][d];
	}
	return 0;
}

/* read feature vectors from a text file */
int readFeatures(char* fname, TypeFeatureList *features) {
    FILE *f;
    long d, i, vectDim, nbFrame;
    double x;
    
    /* Verify if features filename is already for loading */
    
    if ((f = fopen(fname, "r")) == NULL) 
    {
	fprintf(stderr, "ERROR: Features file: %s cannot access! Please verify it! \n", fname);
	exit(EXIT_FAILURE);
    }
    
    if (fscanf(f, "%ld %ld", &nbFrame, &vectDim) != 2)
    {
    	perror("fscanf");
	return 1;
    }

    if (isVerbose) printf("Total frames: %ld   Feature vector size: %ld\n", nbFrame, vectDim);

    createFeatures(features, vectDim, nbFrame);
    
    /* on lit tout le fichier de donnees */
    for(i=0; i<features->nbFrame; i++) {
  	for (d=0; d<features->vectDim; d++) {
		if (fscanf(f, "%lf", &x) != 1) {
		    perror("fscanf");
		    return 1;
		}
		features->featureVectors[i][d] = x;
	}
    }
    
    if (fclose(f)) {
	perror("fclose");
	return 1;
    }
    return 0;
}

/* print feature vectors to the screen */
void printFeatures(TypeFeatureList *features) {
	long d, i;

	printf("Total frames: %ld   Feature vector size: %ld\n", features->nbFrame, features->vectDim);	
	for(i=0; i<features->nbFrame; i++) {
	        for (d=0; d<features->vectDim; d++) 
			printf("%lf\t", features->featureVectors[i][d]);
		printf("\n");
	}

}

/* centering of feature vectors */
void centeringFeatures(TypeFeatureList *features) {
	long d, i;
	double mean;
	for (d=0; d<features->vectDim; d++) {
		mean = 0.0;
		for(i=0; i<features->nbFrame; i++) {
			mean += features->featureVectors[i][d];
		}		
		mean /= features->nbFrame;
		for(i=0; i<features->nbFrame; i++) {
			features->featureVectors[i][d] -= mean;
		}
	}
}

/* decaling of a feature vectors */
void decalingFeatures(TypeFeatureList *features, long step) {
	double *tmp;
	long i, d, abs_step;

	abs_step = abs(step);
	tmp = (double *) malloc (sizeof(double) * (abs_step));
	if (tmp == NULL) {
		fprintf(stderr, "ERROR: Couldn't allocate memory for data buffer (tmp)\n");
		exit(EXIT_FAILURE);
	}

	//Left decaling
	if (step < 0) {
		for (d=0; d<features->vectDim; d++) {
			for (i=0; i<abs_step; i++) 
				tmp[i] = features->featureVectors[i][d];
			for (i=0; i<features->nbFrame-abs_step; i++)
				features->featureVectors[i][d] = features->featureVectors[i+abs_step][d];
			for (i=0; i<abs_step; i++)
				features->featureVectors[features->nbFrame-abs_step+i][d] = tmp[i];
		}
	} else { //Right decaling
		for (d=0; d<features->vectDim; d++) {
			for (i=0; i<abs_step; i++) 
				tmp[i] = features->featureVectors[features->nbFrame-abs_step+i][d];
			for (i=features->nbFrame-1; i>=abs_step; i--)
				features->featureVectors[i][d] = features->featureVectors[i-abs_step][d];
			for (i=0; i<abs_step; i++)
				features->featureVectors[i][d] = tmp[i];
		}
	}
	
	free(tmp);
}

/* Co-inertia Analysis (COIA) of two series of feature vectors
      CoI = trace(X*diag(Dp)*X'*diag(Dn)*Y*diag(Dq)*Y'*diag(Dn));
      for details, see : Dray, Chessel, Thioulouse, "Co-inertia analysis and the linking of ecological data tables" */
double calculateCoInertia(TypeFeatureList *X0, TypeFeatureList *Y0, TypeFeatureList *Weight) {
	TypeFeatureList *X, *Y, tmpMat;
	double CoI, x;
	long i, j, d;

	//Create a temporal square matrix
	initFeatures(&tmpMat);
	createFeatures(&tmpMat, X0->nbFrame, X0->nbFrame);

	if (X0->vectDim > Y0->vectDim) {
		X = Y0;
		Y = X0;
	} else {
		X = X0;
		Y = Y0;
	}

	//tmpMat = X*diag(X->vectDim)
	for (d=0; d<X->vectDim; d++) {
		for(i=0; i<X->nbFrame; i++) tmpMat.featureVectors[i][d] = X->featureVectors[i][d]/X->vectDim;
	}

	//tmpMat = tmpMat*X'
	/*for(j=0; j<X->nbFrame; j++) {
		for(i=0; i<tmpMat.nbFrame; i++) {
			x = 0.0;
			for (d=0; d<X->vectDim; d++) {
				x += X[i
			}
		}
	}*/
	
	
	return CoI;
}

/*
void initMatix(TypeMatix *matrix) {
	matrix->nbRow = matrix->nbCol = 0;
	matrix->mat = NULL;
}

int createMatrix(TypeMatrix *matrix, long nbRow, long nbCol) {
        long i;
	
	matrix->nbRow = nbRow;
	matrix->nbCol = nbCol;
	matrix->mat = (double **) malloc (sizeof(double *) * matrix->nbRow);
	if (matrix->mat == NULL) {
		fprintf(stderr, "ERROR: Couldn't allocate memory for data buffer ((matrix->mat)\n");
		exit(EXIT_FAILURE);
	}
	for(i=0; i<matrix->nbRow; i++) {
		(matrix->mat[i] = (double *) malloc (sizeof(double) *  matrix->nbCo);
		if (matrix->mat[i] == NULL) {
			fprintf(stderr, "ERROR: Couldn't allocate memory for data buffer (matrix->mat[%ld])\n", i);
			exit(EXIT_FAILURE);
		}
	}
	return 0;
}

void deleteMatix(TypeMatrix *matrix) {
	 matrix->nbRow = matrix->nbCol = 0;
	if (matrix->mat != NULL) {
		free(matrix->mat);
		matrix->mat = NULL;
	}
}

int copyMatix(TypeMatix *source_mat, TypeMatix *target_mat) {
        long i, j;
	if ((target_mat->mat != NULL) && (target_mat->nbRow != 0))
		deleteMatix(target_mat);
		createMatrix(target_mat, source_mat->nbRow, source_mat->nbCol);
	for(i=0; i<source_mat->nbRow; i++) {
		for (j=0; j<source_mat->nbCol; j++)
			target_mat->mat[i][j] = source_mat->mat[i][j];
	}
	return 0;
}
*/







/*void multipleWeight(TypeFeatureList *Xs, TypeFeatureList *Vect, TypeFeatureList *Xt) {
	long i, d;
	for (
}*/

/* calculation of the determinant of a covariance matrix */
double getDeterminant(double *cov, long dim, int gaussType) {
    double det, detL;
    long d;
    double sum;
    long i, j, k;
    double *L;
    
    det = detL = 1.0;
    
    if (gaussType == DIAG) {
	for (d=0; d<dim; d++) {
	    det *= cov[d];
	}
    }
    else {  //FULL
	/* Compute cholesky decomposition of a symetric matrix: Cov = L*L' */
	/* det(cov) = det(L)*det(L') = det(L)*det(L) */
	/* Theory: http://www.math-linux.com/spip.php?article35 */
	/* Code: Numerical Recipes Online Book, Section 2.9 */
	
	L = (double *) malloc (sizeof(double) * dim * dim);
	if (L == NULL) {
	    fprintf(stderr, "ERROR: Couldn't allocate memory for data buffer (L)\n");
	    exit(EXIT_FAILURE);
	}
	
	for (i=0; i<dim; i++) {
	    for (j=i; j<dim; j++) {
		L[i*dim+j] = cov[i*dim+j];
		if (i != j) L[j*dim+i] = L[i*dim+j]; //SYM: L(i, j) = L(j, i);
	    }
	}
	
	/*for (i=0; i<dim; i++) {
			for (j=0; j<dim; j++) {
				printf("%8.2f", L[i*dim+j]);
			}
			printf("\n");
		}
		exit(0);*/
	
	for (i=0; i<dim; i++) {
	    for (j=i; j<dim; j++) {
		sum = L[i*dim + j];
		for (k=i-1; k>=0; k--) sum -= L[i*dim+k] * L[j*dim+k];
		if (i == j) {
		    if (sum <= 0.0) { 
			fprintf(stderr, "WARNING: Cholesky decomposition failed in getDeterminant function!\n");
		    }
		    L[i*dim+i] = sqrt(sum); // L[i, i] = sqrt(sum)
		    detL *= L[i*dim+i];
		} else L[j*dim+i] = sum / L[i*dim+i]; //L[j, i] = sum / L[i, i]
	    }
	}
	
	free(L);
	
	det = detL * detL;
    }
    
    return det;
}

/* inversion of a square matrix by using Gauss-Jordan elimination method */
int getInverse(double *m, double *mInv, long dim, bool isSymmetric) {
    double *tmp, akk, aik;
    long i, k, j;
    
    tmp = (double *) malloc (sizeof(double) * dim * dim);
    
    for (i=0; i<dim; i++) {
	for (j=0; j<dim; j++) {
	    if ((j<i) && isSymmetric) tmp[i*dim+j] = m[j*dim+i]; //Symetrical Matrix
	    else tmp[i*dim+j] = m[i*dim+j];
	    
	    if (i==j) mInv[i*dim+j] = 1.0;
	    else  mInv[i*dim+j] = 0.0;
	}
    }
    
    for (k=0; k<dim; k++) {
	akk = tmp[k*dim+k];
	if (akk == 0.0) { 
		fprintf(stderr, "ERROR: Matrix is not inversible!\n");
		exit(EXIT_FAILURE);
	}
	for (j=0; j<dim; j++) {
	    tmp[k*dim+j] = (1.0/akk) * tmp[k*dim+j];
	    mInv[k*dim+j] = (1.0/akk) * mInv[k*dim+j];
	}
	
	for (i=0; i<dim; i++) {
	    if (i==k) continue;
	    aik = tmp[i*dim+k];
	    for (j=0; j<dim; j++) {
		tmp[i*dim+j] = tmp[i*dim+j] - aik * tmp[k*dim+j];
		mInv[i*dim+j] = mInv[i*dim+j] - aik * mInv[k*dim+j];
	    }
	}
	
    }
    
    free(tmp);
    return 0;
}

