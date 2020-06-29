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
 * Authors: Marc Sosnick                                                       *
 * Contributors: Kai J. Kohlhoff, William Hsu                                  *
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
 * \file defaults.h
 * \brief Header for storage and retreival of program parameters.
 *
 * Handles program parameters, including defaults and command-line
 * 
 * \author Author: Marc Sosnick Contributors: Kai J. Kohlhoff, William Hsu
 * \date 5/26/10
 * \version 1.0
 */
#ifndef __DEFAULTS__H__
#define __DEFAULTS__H__

#include <iostream>
#include <string>
#include <sstream>
#include "tokens.h"

//using namespace::std;

#define DEFAULTS_DEBUG 1            /**< set to 1 for extra debugging info at runtime*/

#define CAMPAIGN_DEFAULTS_SUCCESS 0         /**< Return value upon success */
#define CAMPAIGN_DEFAULTS_FAILURE 1     /**< Return value upon error. */

#define DEFAULTS_DEFAULT_DETECTDEVICE true  /**< Default value for device detection */
#define DEFAULTS_DEFAULT_TIMEROUTPUT true   /**< Default value for outputting timer data */
#define DEFAULTS_DEFAULT_DATAPOINTS (0)     /**< Default number of datapoints */
#define DEFAULTS_DEFAULT_DIMENSIONS (0)     /**< Default number of dimensions. */
#define DEFAULTS_DEFAULT_CLUSTERS (0)       /**< Default number of clusters. */
#define DEFAULTS_DEFAULT_INPUTFILENAME ""   /**< Default input filename. */
#define DEFAULTS_DEFAULT_DEVICE (0)     /**< Default device. */
#define DEFAULTS_DEFAULT_LISTDEVICES false  /**< Default value for listing devices to console */

/** 
 * \brief Storage and retreival of program and algorithm parameters.
 * 
 * Defaults stores program and algorithm defaults for the currently
 * implemented algorithms.  Interprets parsed command line into
 * program settings.
 * 
 * \author Author: Marc Sosnick, Contributors: Kai J. Kohlhoff, William Hsu
 * \date 5/26/10
 * \version 1.0
 */
class Defaults{
    
    Tokens commandLineParameters;   /**< Command line parameters */
    
    bool    detectDevice;       /**< If true detects CUDA-enable devices */
    bool    timerOutput;        /**< If true, outputs timer to stdio */
    int dataPoints;     /**< Default number of datapoints */
    int dimensions;     /**< Default number of dimensions */
    int clusters;       /**< Default number of clusters */
    std::string  inputFileName;      /**< Default input filename */
    int device;         /**< Current device */
    bool    listDevices;        /**< If true, lists devices to console and exit */
    
public:
    Defaults();
    ~Defaults();
    
    /**
     * \brief Constructor with command line parsing.
     * constructor initializes object with defaults, parses command line and changes
     * default values appropriately.
     */
    Defaults(int argc, const char **argv);
    
    /**
     * \brief Constructor with command line parsing and clustering algorithm output.
     * Constructor initializes object with defaults, parses command line and changes
     * default values appropriately.  Allows passing of a code that identifies which
     * clustering algorithm is being used.  This can be used later for
     * customizing parameters from tokenizer.
     * \param argc Count of command line arguments.
     * \param argv Array of pointers to command line arguments.
     * \param algo Code for algorithm type.
     */
    Defaults(int argc, const char **argv, std::string algo);
    
    
    /**
     *  \brief Outputs usage to stdio.
     *  \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *          CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int     printUsage();
    
    
    /**
     *  \brief Outputs error message and usage to stdio, exiting with error 1.
     *  \param errorMessage Error message to display
     */
    void    printUsageAndExit(std::string errorMessage);
    
    
    /**
     * \brief Reports if input coming from stdio.
     * \return true if input from stdio, false otherwise.
     */
    bool    isStdIoInput();
    
    
    /**
     * \brief Reports if timer data will be output to stdio.
     * \return true if data will be output, false otherwise.
     */ 
    bool    getTimerOutput();
    
    /**
     * \brief Sets if timer data will be output to stdio.
     * \param valueIn True to output data, false to supress data output.
     * \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *         CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int setTimerOutput(bool valueIn);
    
    
    /**
     * \brief Gets the number of data points.
     * \return Number of datapoints in current dataset.
     */
    int     getDataPoints();
    
    /**
     * \brief Sets the number of data points.
     * \param valueIn number of datapoints.
     * \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *         CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int setDataPoints(int valueIn);
    
    /**
     * \brief Gets the number of dimensions.
     * \return number of dimensions
     */ 
    int     getDimensions();
    
    /**
     * \brief Sets the number of dimensions.
     * \param valueIn number of dimensions
     * \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *         CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int     setDimensions(int valueIn);
    
    
    /**
     * \brief Gets the number of clusters.
     * \return number of clusters
     */ 
    int     getClusters();
    
    /**
     * \brief  Sets the number of clusters.
     * \param  valueIn Number of clusters.
     * \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *         CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int     setClusters(int valueIn);
    
    
    /**
     * \brief Gets the data file name.
     * \return the name of the data file
     */ 
    std::string  getInputFileName();
    
    /**
     * \brief Sets the input data file name
     * \param valueIn input data file name
     * \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *         CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int     setInputFileName(std::string valueIn);
    
    
    /**
     * \brief Gets detect device
     * If detect device is true, the system will be searched
     * for the fastest device.
     * \return State of detect device
     */ 
    bool    getDetectDevice();
    
    /**
     * \brief Sets detect device
     * \param valueIn true sets device detection on, false turns off device detection
     * \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *         CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int     setDetectDevice(bool valueIn);
    
    
    /**
     * \brief Gets list devices
     * If list devices is true, a list of the detected CUDA-enabled devices
     * will be printed to sdtdout, and the program will exit.
     * \return list devices
     */ 
    bool    getListDevices();
    
    /**
     * \brief Sets list devices.
     * \param valueIn If true, the list will be output to stdout, if false no list will be output.
     * \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *         CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int     setListDevices(bool valueIn);
    
    
    /**
     * \brief Gets the current CUDA-enabled device.
     * \return the current CUDA-enabled device
     */ 
    int     getDevice();
    
    /**
     * \brief Sets the current CUDA-enabled device.
     * \param valueIn the CUDA-enabled device address.
     * \return CAMPAIGN_DEFAULTS_SUCCESS upon success, 
     *         CAMPAIGN_DEFAULTS_FAIULRE upon failure.
     */
    int     setDevice(int valueIn);
    
private:
    /**
     * \brief Initializes object
     */
    int     init();
    
    /**
     * \brief Parses command line parameters and sets appropriate parameters.
     * 
     * Processes an already parsed command line, looking for and setting appropriate
     * program parameters.  If there is an error in the parameters, the error
     * will be displayed to console, the usage will be displayed, and the program 
     * will exit with error code 
     */
    int processCommandLineParameters();
    
    /**
     * \brief Checks to see that parameter param has no agruments.
     *
     * If parameter does have arguments, prints an error message, usage,
     * and exits with error.
     * \param param parameter to check for arguments
     */
    void    checkHasNoArguments(std::string param);
    
    /**
     * \brief Checks to see if parameter param has agruments.
     *
     * If parameter doesn't have arguments, prints an error message, usage,
     * and exits with error.
     * \param param parameter to check for arguments
     */
    void    checkHasArguments(std::string param);
    
    /** 
     * \brief Checks to see if parameter param has an integer argument.
     * 
     * If parameter doesn't have an argument, exits with error.
     * If parameter does have an argument, but it is not an integer, 
     * exits with error.
     * \param param parameter to check
     */
    void    checkArgumentHasInteger(std::string param);
};



#endif
