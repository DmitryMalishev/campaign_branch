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
  * \file kcentersGPU.h
  * \brief A CUDA K-centers implementation
  *
  * Implements parallel K-centers clustering on the GPU
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
#define THREADSPERBLOCK 8 /** < Threads per block (tpb) */
#define FLOAT_TYPE float         /** < Precision of floating point numbers */

#undef _GLIBCXX_ATOMIC_BUILTINS

#include <CL/sycl.hpp>
#include <dpct/dpct.hpp>
#include <iostream>
#include <cfloat>
#include "dataio.h"
#include "defaults.h"
#endif


//using namespace std;

/**
 * \brief Parallel algorithm, finds maximum value in first array and key for that value in second array
 * Runtime O(log(BLOCKSIZE)) = O(1)
 * Works for up to 1024 elements in arrays
 * Called from within a kernel, will be inlined
 *
 * \param tid Thread ID
 * \param s_A Array in shared memory
 * \param s_B Array in shared memory
 * \return Maximum value at first position of s_A and corresponding key at first postion of s_B
 */
template <unsigned int BLOCKSIZE, class T, class U>
/* DPCT_ORIG __device__ static void parallelMax(int tid, T *s_A, U *s_B);*/
static void parallelMax(int tid, T* s_A, U* s_B, sycl::nd_item<3> item_ct1);

/**
 * \brief Checks if data points have to be reassigned to current centroid
 * Runtime O(D*N)
 * Each block of threads processes THREADSPERBLOCK data points
 *
 * \param N Number of data points
 * \param D Number of dimensions
 * \param iter Current iteration
 * \param centroid Index of current centroid
 * \param X Array of data point positions
 * \param CTR Coordinates of current centroid
 * \param DIST Distances of data points to their assigned centroids
 * \param ASSIGN Assignments of data points to clusters
 * \param MAXDIST Maximum distance of a data point to its centroid for each block
 * \param MAXID Index of data point furthes away  from its centroid for each block
 * \return Updated values of DIST, ASSIGN, MAXDIST, and MAXID
 **/
 /* DPCT_ORIG __global__ void checkCentroid_CUDA(int N, int D, int iter, int
  * centroid, FLOAT_TYPE *X, FLOAT_TYPE *CTR, FLOAT_TYPE *DIST, int *ASSIGN,
  * FLOAT_TYPE *MAXDIST, int *MAXID);*/
void checkCentroid_sycl(int N, int D, int iter, int centroid, FLOAT_TYPE* X,
    FLOAT_TYPE* CTR, FLOAT_TYPE* DIST, int* ASSIGN,
    FLOAT_TYPE* MAXDIST, int* MAXID,
    sycl::nd_item<3> item_ct1, uint8_t* dpct_local);

/**
 * \brief Runs k-centers on the GPU. Requires CUDA-enabled graphics processor
 * Note: distances should all be FLT_MAX the first time this method is called
 * Runtime O(D*K*N)
 *
 * \param N Number of data points
 * \param K Number of clusters
 * \param D Number of dimensions
 * \param x Clustering input data
 * \param assign Initial assignments of data points to clusters
 * \param dist Distances from each data point to centroid, initialize with FLT_MAX
 * \param centroids Indices of centroids
 * \param seed Index of first centroid
 * \return Updated values of assign, dist, and centroid indices in centroids
 */
void kcentersGPU(int N, int K, int D, FLOAT_TYPE* x, int* assign, FLOAT_TYPE* dist, int* centroids, int seed, DataIO* data);




