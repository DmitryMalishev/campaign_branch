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
  * \file kcentersCPU.cc
  * \brief Test/example file for k-centers clustering on the CPU
  *
  * Test/example file for k-centers clustering on the CPU.
  *
  * \author Author: Kai J. Kohlhoff, Contributors: Marc Sosnick, William Hsu
  * \date 12/2/2010
  * \version 1.0
  **/

#include "defaults.h"
#include "dataio.h"
#include "kcentersCPU.h"
#include "kcentersGPU_Intel.h"
#include "kcentersGPU_Nvidia.h"

//using namespace std;

  /**
   * \brief Main for testing
   */
int main(int argc, const char* argv[])
{
    Defaults* defaults = new Defaults(argc, argv, "kc");

    const int seed = 0;
    DataIO* data = new DataIO;

    float* x = data->readData(defaults->getInputFileName().c_str());
    int N = data->getNumElements();
    int K = data->getNumClusters();
    int D = data->getDimensions();
    FLOAT_TYPE* dist = (FLOAT_TYPE*)malloc(sizeof(FLOAT_TYPE) * N); // cluster distances
    for (int i = 0; i < N; i++) dist[i] = FLT_MAX;
    int* centroids = (int*)malloc(sizeof(int) * K);  // centroid indices
    memset(centroids, 0, sizeof(int) * K);
    int* assign = (int*)malloc(sizeof(int) * N);     // assignments
    memset(assign, seed, sizeof(int) * N);

    // --> do clustering on CPU 
    //kcentersCPU(N, K, D, x, assign, dist, centroids, seed);

    // --> do clustering on GPU
    kcentersGPU(N, K, D, x, assign, dist, centroids, seed, data);

    // print results
    FLOAT_TYPE* ctr = (FLOAT_TYPE*)malloc(sizeof(FLOAT_TYPE) * K * D); // cluster centers
    // for each centroid
    for (int i = 0; i < K; i++)
        // for each dimension
        for (int d = 0; d < D; d++)
            // collect centroid coordinates
            ctr[i * D + d] = x[centroids[i] * D + d];
    data->printClusters(N, K, D, x, ctr, assign);
    free(ctr); ctr = NULL;

    // free memory
    free(x);
    free(dist);
    free(ctr);
    free(assign);

    // done
    std::cout << "Done clustering" << std::endl;

    return 0;
}