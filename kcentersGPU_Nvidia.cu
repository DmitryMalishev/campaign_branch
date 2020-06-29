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
 * Authors:  Kai J. Kohlhoff                                                   *
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
 * \file kcentersGPU.cu
 * \brief A CUDA K-centers implementation
 *
 * Implements parallel K-centers clustering on the GPU
 * 
 * \author Author: Kai J. Kohlhoff, Contributors: Marc Sosnick, William Hsu
 * \date 12/2/2010
 * \version 1.0
 **/

#include "kcentersGPU.h"

using namespace std;


template <unsigned int BLOCKSIZE, class T, class U>
__device__ static void parallelMax(int tid, T *s_A, U *s_B)
{
    if (BLOCKSIZE >= 1024) { if (tid < 512 && s_A[tid + 512] > s_A[tid]) { s_A[tid] = s_A[tid + 512]; s_B[tid] = s_B[tid + 512]; } __syncthreads(); }
    if (BLOCKSIZE >=  512) { if (tid < 256 && s_A[tid + 256] > s_A[tid]) { s_A[tid] = s_A[tid + 256]; s_B[tid] = s_B[tid + 256]; } __syncthreads(); }
    if (BLOCKSIZE >=  256) { if (tid < 128 && s_A[tid + 128] > s_A[tid]) { s_A[tid] = s_A[tid + 128]; s_B[tid] = s_B[tid + 128]; } __syncthreads(); }
    if (BLOCKSIZE >=  128) { if (tid <  64 && s_A[tid +  64] > s_A[tid]) { s_A[tid] = s_A[tid +  64]; s_B[tid] = s_B[tid +  64]; } __syncthreads(); }
    
    if (tid < 32)
    {
        if (BLOCKSIZE >= 64) { if (s_A[tid + 32] > s_A[tid]) { s_A[tid] = s_A[tid + 32]; s_B[tid] = s_B[tid + 32]; } }
        if (BLOCKSIZE >= 32) { if (s_A[tid + 16] > s_A[tid]) { s_A[tid] = s_A[tid + 16]; s_B[tid] = s_B[tid + 16]; } }
        if (BLOCKSIZE >= 16) { if (s_A[tid +  8] > s_A[tid]) { s_A[tid] = s_A[tid +  8]; s_B[tid] = s_B[tid +  8]; } }
        if (BLOCKSIZE >=  8) { if (s_A[tid +  4] > s_A[tid]) { s_A[tid] = s_A[tid +  4]; s_B[tid] = s_B[tid +  4]; } }
        if (BLOCKSIZE >=  4) { if (s_A[tid +  2] > s_A[tid]) { s_A[tid] = s_A[tid +  2]; s_B[tid] = s_B[tid +  2]; } }
        if (BLOCKSIZE >=  2) { if (s_A[tid +  1] > s_A[tid]) { s_A[tid] = s_A[tid +  1]; s_B[tid] = s_B[tid +  1]; } }
    }
}

/**
 * \brief Calculates distance contribution along one dimension
 * It is often more efficient to calculate distances component-wise 
 * (i.e. dimension by dimension) to make optimal use of fast memory
 * (e.g. coalesced memory reads, use of shared memory)
 * Using component-wise calculation requires, for some metrics, the
 * call to distanceFinalize() to complete the distance calculation
 * \param elementA Coordinate A
 * \param elementB Coordinate B
 * \return 1D component of distance between element A and element B
 */
template <class T>
__device__ static T distanceComponentGPU(T *elementA, T *elementB)
{
    T dist = 0.0f;
#ifdef CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED
    // Euclidean distance squared
    dist = elementA[0] - elementB[0];
    dist = dist * dist;
#elif defined(CAMPAIGN_DISTANCE_EUCLIDEAN) // slow because of square root, use of CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED preferred
    // Euclidean distance
    dist = elementA[0] - elementB[0];
    dist = dist * dist;
#elif defined(CAMPAIGN_DISTANCE_MANHATTAN)
    // Manhattan distance
    dist = fabs(elementA[0] - elementB[0]); 
#elif defined(CAMPAIGN_DISTANCE_CHEBYSHEV)
    // Chebyshev distance
    dist = fabs(elementA[0] - elementB[0]);
#else
#error "No distance defined, add #define x to program's header with x = CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED, CAMPAIGN_DISTANCE_EUCLIDEAN, CAMPAIGN_DISTANCE_MANHATTAN, or CAMPAIGN_DISTANCE_CHEBYSHEV");
#endif
    return dist; 
}

/**
 * \brief Finalizes component-wise distance calculation
 * \param D Number of dimensions
 * \param components Pre-computed distance component(s)
 * \return Final distance
 */
template <class T>
__device__ static T distanceFinalizeGPU(int D, T *components)
{
    T dist = 0.0f;
#ifdef CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED
    // Euclidean distance squared
    for (unsigned int cnt = 0; cnt < D; cnt++) dist += components[cnt];
#elif defined(CAMPAIGN_DISTANCE_EUCLIDEAN) // slow because of square root, use of CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED preferred
    // Euclidean distance
    for (unsigned int cnt = 0; cnt < D; cnt++) dist += components[cnt];
    dist = sqrt(dist);
#elif defined(CAMPAIGN_DISTANCE_MANHATTAN)
    // Manhattan distance
    for (unsigned int cnt = 0; cnt < D; cnt++) dist += components[cnt]; 
#elif defined(CAMPAIGN_DISTANCE_CHEBYSHEV)
    // Chebyshev distance
    for (unsigned int cnt = 0; cnt < D; cnt++) dist = max(dist, components[cnt]);
#else
#error "No distance defined, add #define x to program's header with x = CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED, CAMPAIGN_DISTANCE_EUCLIDEAN, CAMPAIGN_DISTANCE_MANHATTAN, or CAMPAIGN_DISTANCE_CHEBYSHEV");
#endif
    return dist;
}

/**
 * \brief Computes distance between element A and element B
 * For best performance elements A and B should be in shared memory
 * \param D Number of dimensions
 * \param elementA Position data of element A
 * \param elementB Position data of element B 
 * \return Distance between element A and element B
 */
template <class T>
__device__ static T distanceGPU(int D, T *elementA, T *elementB) 
{
    T dist = 0.0f;
#ifdef CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED
    // Euclidean distance squared
    for (unsigned int cnt = 0; cnt < D; cnt++)
    {
        T di = (elementA[cnt] - elementB[cnt]);
        dist += di * di;
    }
#elif defined(CAMPAIGN_DISTANCE_EUCLIDEAN) // slow because of square root, use of CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED preferred
    // Euclidean distance
    for (unsigned int cnt = 0; cnt < D; cnt++)
    {
        T di = (elementA[cnt] - elementB[cnt]);
        dist += di * di;
    }
    dist = sqrt(dist);
#elif defined(CAMPAIGN_DISTANCE_MANHATTAN)
    // Manhattan distance
    for (unsigned int cnt = 0; cnt < D; cnt++)
    {
        dist += fabs(elementA[cnt] - elementB[cnt]); 
    }
#elif defined(CAMPAIGN_DISTANCE_CHEBYSHEV)
    // Chebyshev distance
    for (unsigned int cnt = 0; cnt < D; cnt++)
    {
        dist = max(dist, fabs(elementA[cnt] - elementB[cnt]));
    }
#else
#error "No distance defined, add #define x to program's header with x = CAMPAIGN_DISTANCE_EUCLIDEAN_SQUARED, CAMPAIGN_DISTANCE_EUCLIDEAN, CAMPAIGN_DISTANCE_MANHATTAN, or CAMPAIGN_DISTANCE_CHEBYSHEV");
#endif
    return dist; 
    
}

__global__ void checkCentroid_CUDA(int N, int D, int iter, int centroid, FLOAT_TYPE *X, FLOAT_TYPE *CTR, FLOAT_TYPE *DIST, int *ASSIGN, FLOAT_TYPE *MAXDIST, int *MAXID)
{
    extern __shared__ FLOAT_TYPE array[];                  // shared memory
    FLOAT_TYPE *s_dist    = (FLOAT_TYPE*) array;                // tpb distances
    int   *s_ID      = (int*)   &s_dist[blockDim.x];  // tpb IDs
    FLOAT_TYPE *s_ctr     = (FLOAT_TYPE*) &s_ID[blockDim.x];    // tpb centroid components
    
    unsigned int tid = threadIdx.x;                   // thread ID in block
    unsigned int t   = blockIdx.x * blockDim.x + tid; // global thread ID
    
    s_dist[tid] = 0.0;
    s_ID  [tid] = t;
    if (t < N)
    {
        // compute distance
        FLOAT_TYPE dist    = 0.0;
        int   offsetD = 0;
        // process each centroid component
        while (offsetD < D)
        {
            // read up to tpb components from global to shared memory
            if (offsetD + tid < D) s_ctr[tid] = CTR[offsetD + tid];
            __syncthreads();
            // compute distance in each dimension separately and build sum
            for (int d = 0; d < min(blockDim.x, D - offsetD); d++)
            {
                dist += distanceComponentGPU(s_ctr + d - offsetD, X + (offsetD + d) * N + t);
            }
            offsetD += blockDim.x;
            __syncthreads();
        }
        dist = distanceFinalizeGPU<FLOAT_TYPE>(1, &dist);
        // if new distance smaller then update assignment
        FLOAT_TYPE currDist = DIST[t];
        if (dist < currDist)
        {
            DIST  [t] = currDist = dist;
            ASSIGN[t] = iter;
        }
        s_dist[tid] = currDist;
    }
    __syncthreads();
    // find max distance of data point to centroid and its index in block
    parallelMax<THREADSPERBLOCK>(tid, s_dist, s_ID);
    // write maximum distance and index to global mem
    if (tid == 0) 
    {
        MAXDIST[blockIdx.x] = s_dist[tid];
        MAXID  [blockIdx.x] = s_ID[tid];
    }
}



void kcentersGPU(int N, int K, int D, FLOAT_TYPE *x, int *assign, FLOAT_TYPE *dist, int *centroids, int seed, DataIO *data)
{
    // CUDA kernel parameters
    int numBlock = (int) ceil((FLOAT_TYPE) N / (FLOAT_TYPE) THREADSPERBLOCK);
    dim3 block(THREADSPERBLOCK);
    dim3 grid(numBlock);
    int sMem = (2 * sizeof(FLOAT_TYPE) + sizeof(int)) * THREADSPERBLOCK;
    
    // GPU memory pointers, allocate and initialize device memory
    FLOAT_TYPE *dist_d         = data->allocDeviceMemory<FLOAT_TYPE*>(sizeof(FLOAT_TYPE) * N, dist);
    int   *assign_d       = data->allocDeviceMemory<int*>  (sizeof(int)   * N);
    FLOAT_TYPE *x_d            = data->allocDeviceMemory<FLOAT_TYPE*>(sizeof(FLOAT_TYPE) * N * D, x);
    FLOAT_TYPE *ctr_d          = data->allocDeviceMemory<FLOAT_TYPE*>(sizeof(FLOAT_TYPE) * D);
    FLOAT_TYPE *maxDistBlock_d = data->allocDeviceMemory<FLOAT_TYPE*>(sizeof(FLOAT_TYPE) * numBlock);
    int   *maxIdBlock_d   = data->allocDeviceMemory<int*>  (sizeof(int)   * numBlock);
    
    // Initialize host memory
    FLOAT_TYPE *maxDistBlock = (FLOAT_TYPE*) malloc(sizeof(FLOAT_TYPE) * numBlock);
    int   *maxID        = (int*)   malloc(sizeof(int));
    FLOAT_TYPE *ctr          = (FLOAT_TYPE*) malloc(sizeof(FLOAT_TYPE) * D);
    
    int centroid = seed;
    // for each cluster
    for (int k = 0; k < K; ++k)
    {
        // send centroid coordinates for current iteration to device memory
        for (int d = 0; d < D; d++) ctr[d] = x[d * N + centroid];
        cudaMemcpy(ctr_d, ctr, sizeof(FLOAT_TYPE) * D, cudaMemcpyHostToDevice);      
        centroids[k] = centroid;
        
        // for each data point, check if new centroid closer than previous best and if, reassign
        checkCentroid_CUDA<<<grid, block, sMem>>>(N, D, k, centroid, x_d, ctr_d, dist_d, assign_d, maxDistBlock_d, maxIdBlock_d);
        // CUT_CHECK_ERROR("checkCentroid_CUDA() kernel execution failed");
        
        // get next max distance between data point and centroid
        cudaMemcpy(maxDistBlock, maxDistBlock_d, sizeof(FLOAT_TYPE) * numBlock, cudaMemcpyDeviceToHost);        
        int tempMax = 0;
        for (int i = 1; i < numBlock; i++) 
        {
            if (maxDistBlock[i] > maxDistBlock[tempMax]) tempMax = i;
        }
        cudaMemcpy(maxID, maxIdBlock_d + tempMax, sizeof(int), cudaMemcpyDeviceToHost);
        centroid = maxID[0];
    }
    // copy final assignments back to host
    cudaMemcpy(assign, assign_d, sizeof(int) * N, cudaMemcpyDeviceToHost);
    
    // free up memory
    cudaFree(dist_d);
    cudaFree(assign_d);
    cudaFree(x_d);
    cudaFree(ctr_d);
    cudaFree(maxDistBlock_d);
    cudaFree(maxIdBlock_d);
    free(maxDistBlock);
    free(maxID);
    free(ctr);
}


