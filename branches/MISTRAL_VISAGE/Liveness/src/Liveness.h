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
 
#ifndef _liveness_h_
#define _liveness_h_

#include "CommonTools.h"

#define VERSION "1.0.0"

#define MAX_DECALING_BY_SIDE 		10
#define MAX_NATURAL_DELAYED_FRAME	2
#define MIN_NUMBER_FEATURES		40

void version ( void );
void usage ( void );
void erreurOption ( char *option );
double computeLivenessScore(TypeVector *p);

#endif //#ifndef _liveness_h_
