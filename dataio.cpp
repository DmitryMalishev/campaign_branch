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
 * \file dataio.cpp
 * \brief File and stream I/O routines
 *
 * Implements a number of I/O routines for files and streams.
 * 
 * \author Author: Kai J. Kohlhoff, Contributors: Marc Sosnick, William Hsu
 * \date 2/2/2010
 * \version 1.0
 **/

#include <iostream>
#include <fstream>
#include <cfloat>
#include <string>
#include <sstream>
#include <assert.h>
#include "dataio.h"
#include "defaults.h"


using namespace std;

class Internal 
{
private:
    int N; /** < Number of data points */
    int K; /** < Number of cluster centers */
    int D; /** < Number of dimensions */
    int dataSize; /** < Number of elements in data set */
    bool deviceCheck; /** < Determines if device check should be performed */
    bool printTime; /** < Determines if time of execution should be reported */
    float* data; /** < Storage location for clustering input data */
    const char* fileName; /** < Name of data file */
    const char* execName; /** < Name of executable */
    const char* methodType; /** < Type of method, e.g. type of clustering algorithm */
    
public:
    
    /** 
     * Constructor
     */
    Internal()
    {
        N = K = D = dataSize = 0;
        deviceCheck = true;
        printTime = true;
        data = NULL;
    }
    
    /** 
     * Destructor
     */
    ~ Internal()
    {
        delete data;
    }
    
    /**
     * \brief Get number of data points
     * \return Number of data points
     */
    int getNumElements() { return N; };
    
    /**
     * \brief Get number of clusters
     * \return Number of clusters
     */
    int getNumClusters() { return K; };
    
    /**
     * \brief Get number of dimensions
     * \return Number of dimensions
     */
    int getDimensions() { return D; };
    
    /**
     * \brief Get name of executable
     * \return Name of executable
     */
    const char* getExecName()
    {
        return execName;
    }
    
    /**
     * \brief Get name of data file
     * \return Name of data file
     */
    const char* getFileName()
    {
        return fileName;
    }
    
    /**
     * \brief Get type of method
     * \return Type of method
     */
    const char* getMethodType()
    {
        return methodType;
    }
    
    /**
     * \brief Get number of elements in data file
     * \return Number of elements in data file
     */
    int getDataSize()
    {
        return dataSize;
    }
    
    
    /**
     * \brief Set name of executable
     */
    void setExecName(const char* en)
    {
        execName = en;
    }
    
    /**
     * \brief Set name of data file
     */
    void setFileName(const char* fn)
    {
        fileName = fn;
    }
    
    /**
     * \brief Set type of method
     */
    void setMethodType(const char* mt)
    {
        methodType = mt;
    }
    
    /**
     * \brief Set number of data elements
     */
    void setDataSize(int numData)
    {
        dataSize = numData;
    }
    
    /**
     * \brief Returns pointer to current data, or NULL if none has been read
     * \return An array containing the data
     */
    float* getData()
    {
        return data;
    }  
    
    /**
     * \brief Print parameters
     */
    void printParams()
    {
        std::cout << "Using " << N << " data points, " << K << " clusters, " << D << " dimension(s)" << std::endl;
    }
    
    istream& getline(istream& stream, string& str)
    {
        char ch;
        str.clear();
        while (stream.get(ch) && ch != '\n')
            str.push_back(ch);
        return stream;
    }

    /**
     * \brief Read data from file
     * \param fileName Name of data file
     * \return An array containing the data 
     */
    float* readParsFile(const char* fileName)
    {
        
        string line;
        ifstream infile;
        float pars[3];
        infile.open(fileName, ios::in);
        if (!infile.is_open())
        {
            std::cout << "Error in readParsFile(): Unable to find or open file \"" << fileName << "\"." << std::endl;
            exit(1);
        }
        assert(!infile.fail());
        try
        {
            // read parameters (number of data points, clusters, and dimensions)
            for (int i = 0; i < 3; i++)
            {
                getline(infile, line);
                if (infile.eof()) throw 111;
                istringstream buffer(line.c_str());
                if (!(buffer >> pars[i])) throw 112;
            }
            if (N == 0) N = (int) pars[0]; // if N not otherwise initialized
            if (K == 0) K = (int) pars[1]; // if K not otherwise initialized
            if (D == 0) D = (int) pars[2]; // if D not otherwise initialized
        }
        catch (int e)
        {
            std::cout << "Error in dataIO::readParsFile(): ";
            if (e == 111) std::cout << "reached end of file \"" << fileName << "\" prematurely" << std::endl;
            else if (e == 112) std::cout << "can only read floating point numbers" << std::endl;
            else std::cout << "reading file content failed" << std::endl;
            std::cout << "                             Please check parameters and file format" << std::endl;
            return NULL;
        }
        infile.close();
        assert(!infile.fail());
        return data;
    }
    
    /**
     * \brief Read data from file
     * \param fileName Name of data file
     * \return An array containing the data 
     */
    float* readFile(const char* fileName)
    {
        
        string line;
        ifstream infile;
        float pars[3];
        int numData;
        infile.open(fileName, ios::in);
        if (!infile.is_open())
        {
            std::cout << "Error in readFile(): Unable to find or open file \"" << fileName << "\"." << std::endl;
            exit(1);
        }
        assert(!infile.fail());
        try
        {
            // read parameters (number of data points, clusters, and dimensions)
            for (int i = 0; i < 3; i++)
            {
                getline(infile, line);
                if (infile.eof()) throw 42;
                istringstream buffer(line.c_str());
                if (!(buffer >> pars[i])) throw 1337;
            }
            if (N == 0) N = (int) pars[0]; // if N not otherwise initialized
            if (K == 0) K = (int) pars[1]; // if K not otherwise initialized
            if (D == 0) D = (int) pars[2]; // if D not otherwise initialized
            // read the actual data
            if ((numData = dataSize) == 0) 
            {
                printParams();
                numData = N * D;
            }
            std::cout << "Reading " << numData << " floats" << std::endl;
            data = (float*) malloc(sizeof(float) * numData);
            memset(data, 0, sizeof(float) * numData);
            for (int i = 0; i < numData; i++) 
            {
                getline(infile, line);
                if (infile.eof()) throw 42;
                istringstream buffer(line.c_str());
                if (!(buffer >> data[i])) throw 1337;
            }
        }
        catch (int e)
        {
            std::cout << "Error in dataIO::readFile(): ";
            if (e == 42) std::cout << "reached end of file \"" << fileName << "\" prematurely" << std::endl;
            else if (e == 1337) std::cout << "can only read floating point numbers" << std::endl;
            else std::cout << "reading file content failed" << std::endl;
            std::cout << "                             Please check parameters and file format" << std::endl;
            return NULL;
        }
        infile.close();
        assert(!infile.fail());
        return data;
    }
    
    /**
     * \brief Read data from stream, used internally by readStdIn()
     * \param numbers Number of floats to be read
     * \return An array containing the data
     */
    float* readStdIn(int numbers)
    {
        data = (float*) malloc(sizeof(float) * numbers);
        memset(data, 0, sizeof(float) * numbers);
        for (int i = 0; i < numbers; i++)
        {
            cin >> data[i];
        }
        return data;
    }
    
    /**
     * \brief Read data from stdin
     * \return An array containing the data 
     */
    float* readStdIn()
    {
        int numData = 0;
        std::cout << "No file name provided, reading data from stdin" << std::endl;
        std::cout << "Please input number of data points N:" << std::endl;
        cin >> N;
        std::cout << "Please input number of cluster centers K:" << std::endl;
        cin >> K;
        std::cout << "Please input number of dimensions D:" << std::endl;
        cin >> D;
        if ((numData = dataSize) == 0)
        {
            printParams();
            numData = N * D;
        }
        std::cout << "Reading " << numData << " floats" << std::endl;
        readStdIn(numData);
        return data;
    }
    
}; // End of Internal

/** 
 * Constructor
 */
DataIO::DataIO()
{
    ip = new Internal;
}

/**
 * Destructor
 */
DataIO::~DataIO()
{
    delete ip;
}

/** 
 * \brief Read data from file or stdin. 
 * \param fileName Name of data file, if empty quotes ("") read from stdin
 * \return An array containing the data 
 */
float* DataIO::readData(const char* fileName)
{  
    float* data;
    if (fileName == nullptr) data = ip->readStdIn();
    else data = ip->readFile(fileName);
    return data;
}

/**
 * \brief Returns pointer to current data, or NULL if none has been read
 * \return An array containing the data
 */
float* DataIO::getData()
{
    return ip->getData();
}

/**
 * \brief Get name of data file
 * \return Name of data file
 */
const char* DataIO::getFileName()
{
    return ip->getFileName();
}

/**
 * \brief Get number of data points
 * \return Number of data points
 */
int DataIO::getNumElements() { return ip->getNumElements(); }

/**
 * \brief Get number of clusters
 * \return Number of clusters
 */
int DataIO::getNumClusters() { return ip->getNumClusters(); }

/**
 * \brief Get number of dimensions
 * \return Number of dimensions
 */
int DataIO::getDimensions() { return ip->getDimensions(); }

/**
 * \brief Get number of elements in data file
 * \return Number of elements in data file
 */
int DataIO::getDataSize() { return ip->getDataSize(); }



/**
 * \brief Set number of elements in data file. Overrides default behavio
 * \param numData Number of elements in data file
 */
void DataIO::setDataSize(int numData)
{
    ip->setDataSize(numData);
}

/**
 * \brief Prints list of clusters
 */
void DataIO::printClusters(int numData, int numClust, int numDim, float *data, float *ctr, int *assign)
{
    std::cout << "Data clusters:" << std::endl;
    // Data must have at least one dimension
    if (numDim < 1)
    {
        std::cout << "Error in printClusters(). Dimension D has to be larger than 0" << std::endl;
        return;
    }
    // if data points in one dimension (scalar values)
    else if (numDim == 1)
    {
        for (int k = 0; k < numClust; ++k)
        { 
            std::cout << "Cluster " << k << " (";
            int count = 0;
            for (int n = 0; n < numData; ++n)
            {
                if (assign[n] == k) { std::cout << data[n] << ", "; count++; }
            }
            // remove last two characters to avoid trailing comma
            if (count > 0) std::cout << "\b\b"; 
            if (ctr != NULL) std::cout << ") ctr " << ctr[k] << std::endl;
            else std::cout << ")" << std::endl;
        }
    }
    // if vector length > 1
    else
    {
        // for each cluster center
        for (int i = 0; i < numClust; i++)
        {
            std::cout << "Cluster " << i << " (";
            int count = 0;
            // for each data point
            for (int j = 0; j < numData; j++)
            {
                // if element assigned to current cluster center then print out element
                if (assign[j] == i) 
                {
                    // print out vectors
                    std::cout << "{";
		    // Hsu 9/29/11                    for (int cnt = 0; cnt < numDim; cnt++) std::cout << data[numDim * j + cnt] << ((cnt < numDim-1) ? ", " : "");
                    for (int cnt = 0; cnt < numDim; cnt++) std::cout << data[numData * cnt + j] << ((cnt < numDim-1) ? ", " : ""); // Hsu 9/29/11 new
                    std::cout << "}, ";
                    count++;
                }
            }
            // remove last two characters to avoid trailing comma
            if (count > 0) std::cout << "\b\b";
            if (ctr != NULL) 
            {
                // print cluster center
                std::cout << ") ctr {";
                for (int cnt = 0; cnt < numDim; cnt++) std::cout << ctr[numDim * i + cnt] << ", ";
                std::cout << "\b\b}" << std::endl;;
            }
            else std::cout << ")" << std::endl;
        }
    }
}


/** 
 * Main, intended for unit testing
 */
/*
 int main()
 {
 DataIO* fio = new DataIO;
 fio->readData("../../data/testInput.dat");
 return 0;
 }
 */
