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
 * \file kcentersCPU.cpp
 * \brief A K-centers implementation
 *
 * Implements K-centers clustering
 *
 * \author Author: Kai J. Kohlhoff, Contributors: Marc Sosnick, William Hsu
 * \date 12/2/2010
 * \version 1.0
 **/

#include "kcentersCPU.h"

using namespace std;

// defined globally in campaign.h
// define distance metric, e.g. CAMPAIGN_DISTANCE_MANHATTAN, CAMPAIGN_DISTANCE_EUCLIDEAN_, etc.
#define CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED /** < Type of distance metric */
#define FLOAT_TYPE float         /** < Precision of floating point numbers */

#include "metricsCPU.h"

void kcentersCPU(int N, int K, int D, FLOAT_TYPE *x, int *assign, FLOAT_TYPE *dist, int *centroids, int seed)
{
    int centroid = seed; // initialize first cluster center
    // for each cluster
    for (unsigned int k = 0; k < K; k++)
    {
        FLOAT_TYPE maxDist = -1;
        unsigned int newCtr = 0;
        // store index of current centroid
        centroids[k] = centroid;
        // for each data point
        for (unsigned int n = 0; n < N; n++)
        {
            // get distance
            FLOAT_TYPE distance = 0.0;
            for (unsigned int d = 0; d < D; d++)
            {
                distance += distanceComponentCPU(x + d * N + centroid, x + d * N + n);
            }
            distance = distanceFinalizeCPU<FLOAT_TYPE>(1, &distance);
            // if point is closer to current centroid than previous best 
            if (distance < dist[n])
            {
                assign[n] = k;
                dist[n]   = distance;
            }
            // data point furthest from any centroid will be next centroid
            if (dist[n] > maxDist)
            {
                newCtr  = n;
                maxDist = dist[n];
            }
        }
        centroid = newCtr;
    }
}

