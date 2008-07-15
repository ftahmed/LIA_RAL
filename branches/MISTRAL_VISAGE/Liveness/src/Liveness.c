/***************************************************************************
 *                                                                         *
 *   Copyright (C) 2006 by PAROLE Group - LORIA Laboratory                 *
 *   Address: LORIA - Campus Scientifique - BP 239                         *
 *   Author:  Viet Bac LE                                                  *
 *            54506 Vandoeuvre-l�-Nancy Cedex                             *
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include "Liveness.h"
#include "CommonTools.h"

// --------------------------------------------------------------------------

void version ( void )
{
    printf("Test Liveness - A Video & Audio Synchronisation Test\n");
    printf("version %s. (c) 2008 by LIG Laboratory\n", VERSION);
    printf("contact: viet-bac.le@imag.fr.\n");
    printf("This software is under GPL license.\n");
    printf("-\n");
}

// --------------------------------------------------------------------------

void usage ( void )
{
    version();    
    printf("usage:\n");
    printf("\t-h or --help             	         : print this help.\n");
    printf("\t-v or --version                    : print the version of this speaker segmentation tool.\n");
    printf("\t-a or --audio-feat <feature_file>  : input audio feature file (in TXT format).\n");
    printf("\t-e or --energy <feature_file>      : input energy feature file (in TXT format). These features are extracted from audio file.\n");
    printf("\t-f or --video-feat <feature_file>  : input video feature file (in TXT format).\n");
    printf("\t      --verbose                    : verbose mode (show more intermediate information). default: none.\n");
    printf("---\n");
    
}

// --------------------------------------------------------------------------

void erreurOption ( char *option )
{
    fprintf(stderr, "\nERROR: There is nothing after the option %s\n", option);
    exit ( EXIT_FAILURE );
}

// --------------------------------------------------------------------------
double computeLivenessScore(TypeVector *p) {
	long i;
	double max, min, L = 0.0, pRef, sum, mean = 0.0;
	
	/* Normalisation of p */
	max = min = p->vect[0];
	for (i=1; i<p->vectDim; i++) {
		if (max < p->vect[i]) max = p->vect[i];
		if (min > p->vect[i]) min = p->vect[i];
	}

	if (max == min) {
		fprintf(stderr, "\nERROR: Max and min values of the curve are equal. Is this a flat curve ? \n");
	        exit ( EXIT_FAILURE );
	}

	for (i=0; i<p->vectDim; i++) {
		p->vect[i] = (p->vect[i] - min) / (max - min);
		mean += p->vect[i];
	}
	mean /= p->vectDim;
	//printf("mean = %lf\n", mean);

	/* Search for the maximal value from -MAX_NATURAL_DELAYED_FRAME to 0 */
	pRef = p->vect[MAX_DECALING_BY_SIDE];
	for (i=MAX_DECALING_BY_SIDE-1; i>=MAX_DECALING_BY_SIDE-MAX_NATURAL_DELAYED_FRAME; i--) {
		if (pRef < p->vect[i]) pRef = p->vect[i];
	}

	/* Compute final liveness score */
	sum = 0.0;
	for (i=0; i<p->vectDim; i++) if (p->vect[i] <= pRef) sum += 1.00;
	L = (sum / p->vectDim) * (pRef / mean - 1);

	return L;
}

// --------------------------------------------------------------------------
int main(int argc, char *argv[])
{  
    TypeFeatureList audioFeats, videoFeats, decaleAudioFeats, energyFeats, decaleEnergyFeats;
    unsigned short opt;	
    char audioFeatureFname[MAX_STR], videoFeatureFname[MAX_STR], energyFeatureFname[MAX_STR];
    long i;
    TypeVector CoI;
    double Lscore;
    
    /* INITIAL OPTION PARAMETERS */
    isVerbose = false;
    audioFeatureFname[0] = '\0';
    videoFeatureFname[0] = '\0';

    
    /* EXTRACTION DES OPTIONS */        
    if ( argc < 7 )
    {
	usage();
	return ( EXIT_SUCCESS );
    }
    
    for ( opt=1; opt<argc; opt++ )
    {
	if ( ( strcmp(argv[opt], "-h") == 0 ) || ( strcmp(argv[opt], "--help") == 0 ) )
	{
	    usage();
	    return ( EXIT_SUCCESS );
	}
	else
	    if ( ( strcmp(argv[opt], "-v") == 0 ) || ( strcmp(argv[opt], "--version") == 0 ) )
	    {
	    	version();
	    	return ( EXIT_SUCCESS );
	}
	else
	    if ( strcmp(argv[opt], "--verbose") == 0 )
	    {
		isVerbose = true;
	}
	else
	    if ( ( strcmp(argv[opt], "-a") == 0 ) || ( strcmp(argv[opt], "--audio-feat") == 0 ) )
	    {
	    	opt++;
		if (opt < argc)
			sprintf(audioFeatureFname, argv[opt]);
		else
			erreurOption(argv[--opt]);
	}
	else
	    if ( ( strcmp(argv[opt], "-e") == 0 ) || ( strcmp(argv[opt], "--energy") == 0 ) )
	    {
		opt++;
		if (opt < argc)
			sprintf(energyFeatureFname, argv[opt]);
		else
			erreurOption(argv[--opt]);
	    }
	else
	    if ( ( strcmp(argv[opt], "-f") == 0 ) || ( strcmp(argv[opt], "--video-feat") == 0 ) )
	    {
	    	opt++;
		if (opt < argc)                                
			sprintf(videoFeatureFname, argv[opt]);
		else
			erreurOption(argv[--opt]);
	}
	else {
	    fprintf(stderr, "WARNING: Unknown option %s\n", argv[opt]);
	}
	
    }
    
    version();

    /* Initialisation */
    initFeatures(&audioFeats);
    initFeatures(&energyFeats);
    initFeatures(&videoFeats);
    initFeatures(&decaleAudioFeats);
    initFeatures(&decaleEnergyFeats);
    
    /* Read audio and video features from files */
    readFeatures(audioFeatureFname, &audioFeats);
    readFeatures(energyFeatureFname, &energyFeats);
    readFeatures(videoFeatureFname, &videoFeats);

    /* center the features => mean by dimention is 0 */
    centeringFeatures(&audioFeats);
    centeringFeatures(&videoFeats);
    if (audioFeats.nbFrame != energyFeats.nbFrame) {
    	fprintf(stderr, "\nERROR: Audio features size and energy features size are different\n");
	exit ( EXIT_FAILURE );
    }
    if (audioFeats.nbFrame < MIN_NUMBER_FEATURES) {
    	fprintf(stderr, "\nERROR: Number of audio features is too short! Allowed minimal number of features = %d\n", MIN_NUMBER_FEATURES);
	exit ( EXIT_FAILURE );
    }
    if (videoFeats.nbFrame < MIN_NUMBER_FEATURES) {
    	fprintf(stderr, "\nERROR: Number of video features is too short! Allowed minimal number of features = %d\n", MIN_NUMBER_FEATURES);
	exit ( EXIT_FAILURE );
    }

    if (audioFeats.nbFrame < videoFeats.nbFrame) videoFeats.nbFrame = audioFeats.nbFrame;
    else if (audioFeats.nbFrame > videoFeats.nbFrame) {
    	audioFeats.nbFrame = energyFeats.nbFrame = videoFeats.nbFrame;
    }

    createFeatures(&decaleAudioFeats, audioFeats.vectDim, audioFeats.nbFrame);
    copyFeatures(&audioFeats, &decaleAudioFeats);
    createFeatures(&decaleEnergyFeats, energyFeats.vectDim, energyFeats.nbFrame);
    copyFeatures(&energyFeats, &decaleEnergyFeats);
    //printFeatures(&decaleAudioFeats);
    
    CoI.vectDim = 2*MAX_DECALING_BY_SIDE+1;
    CoI.vect = (double *) malloc (sizeof(double) * CoI.vectDim);
    decalingFeatures(&decaleAudioFeats, -MAX_DECALING_BY_SIDE);
    decalingFeatures(&decaleEnergyFeats, -MAX_DECALING_BY_SIDE);
    for (i=0; i<CoI.vectDim; i++) {    	
	CoI.vect[i] = calculateCoInertia(&videoFeats, &decaleAudioFeats, &decaleEnergyFeats);
	printf("Co-Inertia[%ld] = %lf\n", i-MAX_DECALING_BY_SIDE, CoI.vect[i]);
	decalingFeatures(&decaleAudioFeats, 1);
	decalingFeatures(&decaleEnergyFeats, 1);
    }

    Lscore = computeLivenessScore(&CoI);
    printf("Liveness score by COIA = %lf\n", Lscore);
    
    /* Finalisation */
    deleteFeatures(&audioFeats);
    deleteFeatures(&videoFeats);
    deleteFeatures(&decaleAudioFeats);
    return EXIT_SUCCESS;
}
