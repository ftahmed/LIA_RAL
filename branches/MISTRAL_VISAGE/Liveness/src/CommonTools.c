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

#include "CommonTools.h"

/* init feature vectors */
void initFeatures(TypeFeatureList *features) {
	features->nbFrame = features->vectDim = 0;
	features->nbRow = features->nbRow = 0;
	features->nbRowMax = features->nbRowMax = 0;
	features->featureVectors = NULL;
}

/* create feature vectors and allocate memory */
int createFeatures(TypeFeatureList *features, long vectDim, long nbFrame) {
	long i;

	features->nbFrame = features->nbRow = features->nbRowMax = nbFrame;
	features->vectDim = features->nbCol = features->nbColMax = vectDim;
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
	features->nbRow = features->nbRow = 0;
	features->nbRowMax = features->nbRowMax = 0;
	if (features->featureVectors != NULL) {
		free(features->featureVectors);
		features->featureVectors = NULL;
	}
}

/* copy feature vectors */
int copyFeatures(TypeFeatureList *sourceFeats, TypeFeatureList *targetFeats) {
	long i, d;

	/* if ((targetFeats->featureVectors != NULL) && (targetFeats->nbFrame != 0))
		deleteFeatures(targetFeats);
	createFeatures(targetFeats, sourceFeats->vectDim, sourceFeats->nbFrame); */

	for(i=0; i<sourceFeats->nbFrame; i++) {
	        for (d=0; d<sourceFeats->vectDim; d++) 
			targetFeats->featureVectors[i][d] = sourceFeats->featureVectors[i][d];
	}
	targetFeats->nbFrame = sourceFeats->nbFrame;
	targetFeats->vectDim = sourceFeats->vectDim;
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

/* Update features parameters from matrix parameters */
void updateFeatures(TypeFeatureList *x) {
        x->featureVectors = x->mat;
	x->nbFrame = x->nbRow;
	x->vectDim = x->nbCol;
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
	long i, d, absStep;

	absStep = abs(step);
	tmp = (double *) malloc (sizeof(double) * (absStep));
	if (tmp == NULL) {
		fprintf(stderr, "ERROR: Couldn't allocate memory for data buffer (tmp)\n");
		exit(EXIT_FAILURE);
	}

	//Left decaling
	if (step < 0) {
		for (d=0; d<features->vectDim; d++) {
			for (i=0; i<absStep; i++) 
				tmp[i] = features->featureVectors[i][d];
			for (i=0; i<features->nbFrame-absStep; i++)
				features->featureVectors[i][d] = features->featureVectors[i+absStep][d];
			for (i=0; i<absStep; i++)
				features->featureVectors[features->nbFrame-absStep+i][d] = tmp[i];
		}
	} else { //Right decaling
		for (d=0; d<features->vectDim; d++) {
			for (i=0; i<absStep; i++) 
				tmp[i] = features->featureVectors[features->nbFrame-absStep+i][d];
			for (i=features->nbFrame-1; i>=absStep; i--)
				features->featureVectors[i][d] = features->featureVectors[i-absStep][d];
			for (i=0; i<absStep; i++)
				features->featureVectors[i][d] = tmp[i];
		}
	}
	
	free(tmp);
}

/* Co-inertia Analysis (COIA) of two series of feature vectors
      CoI = trace(X*diag(Dp)*X'*diag(Dn)*Y*diag(Dq)*Y'*diag(Dn));
      for details, see : Dray, Chessel, Thioulouse, "Co-inertia analysis and the linking of ecological data tables" */
double calculateCoInertia(TypeFeatureList *X0, TypeFeatureList *Y0, TypeFeatureList *Weight) {
	TypeFeatureList *X, *Y;
	TypeMatrix tmpMat1, tmpMat2, tmpMat3;
	TypeVector tmpVect;
	double CoI, a, sumE;
	long i;

	if ((X0->nbFrame != Y0->nbFrame) || (X0->nbFrame != Weight->nbFrame)) {
		fprintf(stderr, "ERROR: Sizes of two series of feature vectors X0, Y0 are not unique\n");
                exit(EXIT_FAILURE);
	}

	if (X0->vectDim > Y0->vectDim) {
		X = Y0;
		Y = X0;
	} else {
		X = X0;
		Y = Y0;
	}
	updateMatrix(X);
	updateMatrix(Y);
	
	//Create a temporal square matrix
	initFeatures(&tmpMat1);
	initFeatures(&tmpMat2);
	initFeatures(&tmpMat3);
	createFeatures(&tmpMat1, X0->nbFrame, X0->nbFrame); updateMatrix(&tmpMat1);
	createFeatures(&tmpMat2, X0->nbFrame, X0->nbFrame); updateMatrix(&tmpMat2);
	createFeatures(&tmpMat3, X0->nbFrame, X0->nbFrame); updateMatrix(&tmpMat3);

	//tmpMat1 = X
	copyFeatures(X, &tmpMat1); updateMatrix(&tmpMat1);

	//tmpMat2 = tmpMat1*diag(ones(X->vectDim)/X->vectDim)
	a = 1.0 / X->vectDim;
	valueMultiple(&tmpMat1, a, &tmpMat2);

	//tmpMat1 = tmpMat2*X'
	transposeMatrix(X, &tmpMat3);
	matrixMultiple(&tmpMat2, &tmpMat3, &tmpMat1);

	//tmpMat2 = tmpMat1*diag(Weight)
	tmpVect.vectDim = Weight->nbFrame;
	tmpVect.vect = (double *) malloc (sizeof(double) * tmpVect.vectDim);
	sumE = 0.0;
	for (i=0; i<tmpVect.vectDim; i++) sumE += Weight->featureVectors[i][0];
	for (i=0; i<tmpVect.vectDim; i++) tmpVect.vect[i] = Weight->featureVectors[i][0] / sumE;
	diagMultiple(&tmpMat1, &tmpVect, &tmpMat2);

	//tmpMat1 = tmpMat2*Y;
	matrixMultiple(&tmpMat2, Y, &tmpMat1);

	//tmpMat2 = tmpMat1*diag(ones(Y->vectDim)/Y->vectDim)
	a = 1.0 / Y->vectDim;
	valueMultiple(&tmpMat1, a, &tmpMat2);

	//tmpMat1 = tmpMat2*Y'
	transposeMatrix(Y, &tmpMat3);
	matrixMultiple(&tmpMat2, &tmpMat3, &tmpMat1);

	//tmpMat2 = tmpMat1*diag(Weight)
	diagMultiple(&tmpMat1, &tmpVect, &tmpMat2);

	//Sum of diagonal elements of tmpMat2
	CoI = traceMatrix(&tmpMat2);
	//updateFeatures(&tmpMat2); printFeatures(&tmpMat2);

	deleteFeatures(&tmpMat1);
	deleteFeatures(&tmpMat2);
	deleteFeatures(&tmpMat3);
	free(tmpVect.vect);
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

int copyMatix(TypeMatix *sourceMat, TypeMatix *targetMat) {
        long i, j;
	
	if ((targetMat->mat != NULL) && (targetMat->nbRow != 0))
		deleteMatix(targetMat);
		createMatrix(targetMat, sourceMat->nbRow, sourceMat->nbCol);
	for(i=0; i<sourceMat->nbRow; i++) {
		for (j=0; j<sourceMat->nbCol; j++)
			targetMat->mat[i][j] = sourceMat->mat[i][j];
	}
	return 0;
}
*/

/* Update matrix parameters from features parameters */
void updateMatrix(TypeMatrix *X) {
	X->mat = X->featureVectors;
	X->nbRow = X->nbFrame;
	X->nbCol = X->vectDim;
}

/* Transpose a matrix */
void transposeMatrix(TypeMatrix *X, TypeMatrix *XT) {
	long i, j;

	if ((XT->nbRowMax < X->nbCol) || (XT->nbColMax < X->nbRow)) {
		fprintf(stderr, "ERROR: Sizes of matrix are incompatible in matrix transposing! \n");
                exit(EXIT_FAILURE);
	}
	
	for (i=0; i<X->nbRow; i++) {
	                for (j=0; j<X->nbCol; j++) XT->mat[j][i] = X->mat[i][j];
	}

	XT->nbRow = X->nbCol;
	XT->nbCol = X->nbRow;
}

/* Sum of diagonal elements of a matrix */
double traceMatrix(TypeMatrix *X) {
	double a = 0.0;
	long i;

	if (X->nbRow != X->nbCol) {
		fprintf(stderr, "ERROR: Cannot get the sum of diagonal elements of a non-squared matrix! \n");
		exit(EXIT_FAILURE);
	}

	for (i=0; i<X->nbRow; i++) a += X->mat[i][i];

	return a;
}

/* Multiplying a matrix with another diagonal matrix (generated from a vector) */
void diagMultiple(TypeMatrix *Xs, TypeVector *v, TypeMatrix *Xt) {
	long i, j;

	for (j=0; j<Xs->nbCol; j++) {
		for (i=0; i<Xs->nbRow; i++) Xt->mat[i][j] = Xs->mat[i][j] * v->vect[j];
	}
	Xt->nbRow = Xs->nbRow;
	Xt->nbCol = Xs->nbCol;
}

/* Multiplying a matrix with a number : Xt(m,n) = a * Xs(m,n) */
void valueMultiple(TypeMatrix *Xs, double a, TypeMatrix *Xt) {
	long i, j;

	for (i=0; i<Xs->nbRow; i++) {
		for (j=0; j<Xs->nbCol; j++) {
			Xt->mat[i][j] = Xs->mat[i][j] * a;
		}
	}
	Xt->nbRow = Xs->nbRow;
	Xt->nbCol = Xs->nbCol;
}

/* Multiplying of two matrix : Xt(m,n) = Xs1(m,k) * Xs2(k,n) */
void matrixMultiple(TypeMatrix *Xs1, TypeMatrix *Xs2, TypeMatrix *Xt) {
        long i, j, k;
	double sum;

	if (Xs1->nbCol != Xs2->nbRow) {
		fprintf(stderr, "ERROR: Sizes of matrix are incompatible in matrix multiplying!\n");
		exit(EXIT_FAILURE);
	}
	
	for (i=0; i<Xs1->nbRow; i++) {
		for (j=0; j<Xs2->nbCol; j++) {
			sum = 0.0;
			for (k=0; k<Xs1->nbCol; k++) sum += Xs1->mat[i][k] * Xs2->mat[k][j];
			Xt->mat[i][j] = sum;
		}
	}
	Xt->nbRow = Xs1->nbRow;
	Xt->nbCol = Xs2->nbCol;
}

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

