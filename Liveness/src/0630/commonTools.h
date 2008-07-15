/***************************************************************************
 *   Author:  Viet Bac LE                                                  *
 *   Contact: vietbac.le@loria.fr                                          *
 *                                                                         *
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

#ifndef _commonTools_h_
#define _commonTools_h_

#include <string.h>

#define VERSION "1.0.0"

#define MAX_SEG 5000
#define MAX_RUPTURE 5000
#define MAX_STR 1024

//#define MAX_FLOAT 3.40282347e+38F
#define MAX_DOUBLE 1.79769313486231570e+308
#define MAX_INT 0x7fffffff

#define false 0
#define true 1

#define DIAG 0
#define FULL 1

typedef int bool;

bool isVerbose;

/* data type for feature vectors */
typedef struct TFeatureList {
    double **featureVectors;
    long vectDim, nbFrame;
    double **mat;
    long nbRow, nbCol, nbRowMax, nbColMax;
} TypeFeatureList, TypeMatrix;

typedef struct TVector {
    double *vect;
    long vectDim;
} TypeVector;

/*typedef struct TMatrix {
    double **mat;
    long nbRow, nbCol;
} TypeMatrix;*/


/* init feature vectors */
void initFeatures(TypeFeatureList *features);
/* create feature vectors and allocate memory */
int createFeatures(TypeFeatureList *features, long vectDim, long nbFrame);
/* reset feature vectors and free allocated memory */
void deleteFeatures(TypeFeatureList *features);
/* copy feature vectors */
int copyFeatures(TypeFeatureList *source_feats, TypeFeatureList *target_feats);
/* read feature vectors from a feature file */
int readFeatures(char* fname, TypeFeatureList *features);
/* Update features parameters from matrix parameters */
void updateFeatures(TypeFeatureList *X);
/* print feature vectors to the screen */
void printFeatures(TypeFeatureList *features);
/* centering of feature vectors */
void centeringFeatures(TypeFeatureList *features);
/* decaling of a feature vectors */
void decalingFeatures(TypeFeatureList *features, long step);
/* Co-inertia Analysis (COIA) of two series of feature vectors */
double calculateCoInertia(TypeFeatureList *X0, TypeFeatureList *Y0, TypeFeatureList *Weight);

/* Update matrix parameters from features parameters */
void updateMatrix(TypeMatrix *X);
/* Transpose a matrix */
void transposeMatrix(TypeMatrix *X, TypeMatrix *XT);
/* Sum of diagonal elements of a matrix */
double traceMatrix(TypeMatrix *X);
/* Multiplying a matrix with another diagonal matrix (generated from a vector) */
void diagMultiple(TypeMatrix *Xs, TypeVector *v, TypeMatrix *Xt);
/* Multiplying a matrix with a number : Xt(m,n) = a * Xs(m,n) */
void valueMultiple(TypeMatrix *Xs, double a, TypeMatrix *Xt);
/* Multiplying of two matrix : Xt(m,n) = Xs1(m,k) * Xs2(k,n) */
void matrixMultiple(TypeMatrix *Xs1, TypeMatrix *Xs2, TypeMatrix *Xt);
/* calculation of the determinant of a covariance matrix */
double getDeterminant(double *cov, long dim, int gaussType);
/* inversion of a square matrix by using Gauss-Jordan method */
int getInverse(double *m, double *mInv, long dim, bool isSymmetric);

#endif //#ifndef _commonTools_h_
