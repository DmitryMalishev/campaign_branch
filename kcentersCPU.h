/*
 * --------------------------------------------------------------------------- *
 *                                  CAMPAIGN                                   *
 * --------------------------------------------------------------------------- *
 * This is part of the CAMPAIGN data clustering library originating from       *
 * Simbios, the NIH National Center for Physics-Based Simulation of Biological *
 * Structures at Stanford, funded under the NIH Roadmap for Medical Research,  *
 * grant U54 GM072970 (See https://simtk.org), and the FEATURE Project at      *
 * Stanford, funded under the NIH grant LM05652                                *
 * (See http://feature.stanford.edu/index.php).                                *
 *                                                                             *
 * Portions copyright (c) 2010 Stanford University, Authors, and Contributors. *
 * Authors: Kai J. Kohlhoff                                                    *
 * Contributors: Marc Sosnick, William Hsu                                     *
 *                                                                             *
 * This program is free software: you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                                  *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public        *
 * License for more details.                                                   *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.       *
 * --------------------------------------------------------------------------- *
 */

/* $Id$ */

/**
 * \file kcentersCPU.h
 * \brief A K-centers implementation
 *
 * Implements K-centers clustering
 *
 * \author Author: Kai J. Kohlhoff, Contributors: Marc Sosnick, William Hsu
 * \date 12/2/2010
 * \version 1.0
 **/

#ifdef HAVE_CONFIG_H
#include "../config.h"
#include "../campaign.h"
#else
// defined globally in campaign.h
// define distance metric, e.g. CAMPAIGN_DISTANCE_MANHATTAN, CAMPAIGN_DISTANCE_EUCLIDEAN_, etc.
#define CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED /** < Type of distance metric */
#define FLOAT_TYPE float         /** < Precision of floating point numbers */

#include <iostream>
#include <cfloat>
#include "dataio.h"
#include "defaults.h"
#endif

//using namespace std;

/**
 * \brief Runs K-centers on the CPU
 * Runtime O(D*K*N)
 *
 * \param N Number of data points
 * \param K Number of clusters
 * \param D Number of dimensions
 * \param x Clustering input data
 * \param assign Assignments of data points to clusters
 * \param dist Distances from each data point to centroid, initialize with FLT_MAX
 * \param centroids Centroids
 * \param seed Index of first centroid
 * \return Updated values of assign, dist, and centroid indices in centroids
 */ 
void kcentersCPU(int N, int K, int D, FLOAT_TYPE *x, int *assign, FLOAT_TYPE *dist, int *centroids, int seed);

