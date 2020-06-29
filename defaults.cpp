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
 *  \file defaults.cpp
 *  \brief Storage and retreival of program parameters.
 * 
 *  Handles program parameters, including defaults and command-line parameters
 * 
 *  \author Author: Marc Sosnick, Contributors: Kai J. Kohlhoff, William Hsu
 *  \date 5/26/10
 *  \version 1.0
 */

#include "defaults.h"

using namespace::std;

Defaults::Defaults(){
    if(DEFAULTS_DEBUG) std::cout << "Defaults::Defaults()" << std::endl;
    
    init();
}


Defaults::Defaults(int argc, const char **argv){
    
    if(DEFAULTS_DEBUG) std::cout << "Defaults::Defaults(" << argc << ", " << argv << ")" << std::endl;
    
    init();
    
    commandLineParameters.tokenizeCommandLine(argc, argv);
    
    processCommandLineParameters();
}


Defaults::Defaults(int argc, const char **argv, string algo){
    
    if(DEFAULTS_DEBUG) std::cout << "Defaults::Defaults(" << argc << ", " << argv << ")" << std::endl;
    
    if( algo == "km") std::cout << "K-means" << std::endl;
    else if  ( algo == "kc") std::cout << "K-centers" << std::endl;
    else if (algo == "hier") std::cout << "Hierarchical clustering" << std::endl;
    else if (algo == "som") std::cout << "Self organizing map" << std::endl;
    else if (algo ==  "birch") std::cout << "BIRCH clustering" << std::endl;
    else if (algo ==  "rmsd") std::cout << "Protein backbone rmsd" << std::endl;
    else std::cout<< "Undefined clustering algorithm" << std::endl;
    
    init();
    
    commandLineParameters.tokenizeCommandLine(argc, argv);
    
    processCommandLineParameters();
}

Defaults::~Defaults(){
}



void Defaults::checkHasNoArguments(string param){
    if( commandLineParameters.hasArgument(param)) 
        printUsageAndExit("Error: " + param + " parameter takes no arguments");
}

void Defaults::checkHasArguments(string param){
    if( !commandLineParameters.hasArgument(param)) 
        printUsageAndExit("Error: " + param + " parameter requires an argument");
}

void Defaults::checkArgumentHasInteger(string param){
    if( !commandLineParameters.hasArgument(param))
        printUsageAndExit("Error: " + param + " parameter requires an integer argument");
    
    if( !commandLineParameters.argumentIsNumber(param) )
        printUsageAndExit("Error: unknown argument " + commandLineParameters.getStringArgument(param) + 
                          ". " +  param +  " requires an integer argument."); 
}

int Defaults::processCommandLineParameters(){
    
    if(DEFAULTS_DEBUG) std::cout << "Defaults::processCommandLineParameters()" << std::endl;
    
    // _UNPAIRED parameter means arguments unpaired with paramters. One filename required
    if(commandLineParameters.getArgumentInstances("_UNPAIRED") < 1)
        printUsageAndExit("Error: data file name required.");
    
    // _UNPAIRED parameter means arguments unpaired with paramters. Only one unpaired argument allowed.
    if(commandLineParameters.getArgumentInstances("_UNPAIRED") > 1)
        printUsageAndExit("Error: only one filename allowed on command line.");
    
    // Only one parameter allowed per argument
    for(int i=0; i < commandLineParameters.getParameterCount() ; i++){
        if( commandLineParameters.getArgumentInstances(commandLineParameters.getArgument(i)) > 1)
            printUsageAndExit("Error: duplicate switch: " + commandLineParameters.getArgument(i));
    }
    

    string curParam;
    for (int i = 0; i < commandLineParameters.getParameterCount(); i++) {
        curParam = commandLineParameters.getArgument(i);

        if (curParam == "d") {
            if (commandLineParameters.hasArgument(curParam)) {
                if (commandLineParameters.getStringArgument(curParam) == "true") {
                    setDetectDevice(true);
                }
                else if (commandLineParameters.getStringArgument(curParam) == "false") {
                    setDetectDevice(false);
                }
                else {
                    printUsageAndExit("Erorr: unknown argument for -d: " +
                        commandLineParameters.getStringArgument(curParam));
                }

            }
            else {
                setDetectDevice(true);
            }
        }
        else if (curParam == "detect") {
            checkHasNoArguments(curParam);
            setDetectDevice(true);
        }
        else if (curParam == "nodetect") {
            checkHasNoArguments(curParam);
            if (!commandLineParameters.hasArgument("select")) {
                printUsageAndExit((string)"Error: if you do not detect CUDA device, you must " +
                    "specify which CUDA device to use.  To see list of CUDA devices, " +
                    "use --devices");
            }
            setDetectDevice(false);
        }
        else if (curParam == "devices") {
            checkHasNoArguments(curParam);
            setListDevices(true);
        }
        else if (curParam == "select") {
            checkArgumentHasInteger(curParam);
            setDetectDevice(false);
            setDevice(commandLineParameters.getIntegerArgument(curParam));
        }
        else if (curParam == "t") {
            if (commandLineParameters.hasArgument(curParam)) {
                if (commandLineParameters.getStringArgument(curParam) == "true") {
                    setTimerOutput(true);
                }
                else if (commandLineParameters.getStringArgument(curParam) == "false") {
                    setTimerOutput(false);
                }
                else {
                    printUsageAndExit("Error:  unknown argument for -t: " +
                        commandLineParameters.getStringArgument(curParam));
                }
            }
            else {
                setTimerOutput(true);
            }

        }
        else if (curParam == "timeroutput") {
            checkHasNoArguments(curParam);
            setTimerOutput(true);
        }
        else if (curParam == "notimeroutput") {
            checkHasNoArguments(curParam);
            setTimerOutput(false);
        }
        else if (curParam == "N" || curParam == "dataPoints") {
            checkArgumentHasInteger(curParam);
            setDataPoints(commandLineParameters.getIntegerArgument(curParam));
        }
        else if (curParam == "D" || curParam == "dimensions") {
            checkArgumentHasInteger(curParam);
            setDimensions(commandLineParameters.getIntegerArgument(curParam));
        }
        else if (curParam == "K" || curParam == "clusters") {
            checkArgumentHasInteger(curParam);
            setClusters(commandLineParameters.getIntegerArgument(curParam));
        }
        else if (curParam == "f" || curParam == "file") {
            checkHasArguments(curParam);
            setInputFileName(commandLineParameters.getStringArgument(curParam));
        }
        else if (curParam == "h" || curParam == "help") {
            checkHasNoArguments(curParam);
            printUsage();
            exit(0);
        }
        else if (curParam == "_UNPAIRED") {
            // _UNPAIRED parameter contains input filename 
            setInputFileName(commandLineParameters.getStringArgument(curParam));
        }
        else {
            printUsageAndExit("Error: Unknown command line parameter " + curParam);
        }
    }

    return 0;
}


int Defaults::init(){
    if(DEFAULTS_DEBUG) std::cout << "init()" << std::endl;
    
    detectDevice = DEFAULTS_DEFAULT_DETECTDEVICE;
    timerOutput = DEFAULTS_DEFAULT_TIMEROUTPUT;
    dataPoints = DEFAULTS_DEFAULT_DATAPOINTS;
    dimensions = DEFAULTS_DEFAULT_DIMENSIONS;
    clusters = DEFAULTS_DEFAULT_CLUSTERS;
    inputFileName = DEFAULTS_DEFAULT_INPUTFILENAME;
    device = DEFAULTS_DEFAULT_DEVICE;
    listDevices = DEFAULTS_DEFAULT_LISTDEVICES;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}


void Defaults::printUsageAndExit(string errorMessage){
    std::cout << errorMessage << std::endl;
    printUsage();
    exit(1);
}


int Defaults::printUsage(){
    if(DEFAULTS_DEBUG) std::cout << "printUsage()" << std::endl;
    
    std::cout << "Options: " << std::endl
    << "   -d [bool]            If bool is true, detects CUDA-enabled devices." << std::endl
    << "                        If bool is false, does not detect cuda devices." << std::endl
    << "                        Defaults to true. " << std::endl
    << "                        If CUDA-enabled devices are detected, the fastest" << std::endl
    << "                        device is used.  If CUDA-enabled devices are not " << std::endl
    << "                        detected, device 0 is used. " << std::endl
    << "  --detect              Equivalent to -d true " << std::endl
    << "  --nodetect            Equivalent to -d false " << std::endl;
    std::cout << std::endl
    << "  --devices             Displays a list of CUDA-enabled devices available, " << std::endl
    << "                        along with their index. " << std::endl
    << std::endl
    << "  --select device       Where device is an index of a CUDA-enabled device  " << std::endl
    << "                        displayed using the --devices option. " << std::endl
    << std::endl
    << "   -t [bool]            If bool is true, activates timer output. " << std::endl
    << "                        If bool is false, deactivates timer output. " << std::endl
    << "                        Defaults to true. " << std::endl
    << "  --timeroutput         Equivalent to -t true " << std::endl
    << "  --notimeroutput       Equivalent to -t false " << std::endl
    << std::endl;
    std::cout << "   -N points            Sets number of datapoints to points. Overrides number " << std::endl
    << "                        of points read from file. " << std::endl
    << std::endl 
    << "   -D dims              Sets number of dimensions to dims.  Overrides number of " << std::endl
    << "                        dimensions read from file." << std::endl
    << std::endl
    << "   -K clust             Sets number of clusters to clust.  Overrides number of " << std::endl
    << "                        clusters read from file. " << std::endl
    << std::endl;
    std::cout << "   -f filename " << std::endl
    << "  --file filename       Sets input data file to filename.  Overrides any input " << std::endl
    << "                        from standard input. " << std::endl
    << std::endl;
    std::cout << "   -h " << std::endl
    << "  --help                Displays this help information. " << std::endl;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}

bool Defaults::getDetectDevice(){
    return detectDevice;
}

int Defaults::setDetectDevice(bool valueIn){
    detectDevice = valueIn;
    return CAMPAIGN_DEFAULTS_SUCCESS;
}

bool Defaults::getTimerOutput(){
    return timerOutput;
}

int Defaults::setTimerOutput(bool valueIn){
    timerOutput = valueIn;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}

int Defaults::getDataPoints(){
    return dataPoints;
}

int Defaults::setDataPoints(int valueIn){
    dataPoints = valueIn;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}

int Defaults::getDimensions(){
    return dimensions;
}

int Defaults::setDimensions(int valueIn){
    dimensions = valueIn;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}

int Defaults::getClusters(){
    return clusters;
}

int Defaults::setClusters(int valueIn){
    clusters = valueIn;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}

string Defaults::getInputFileName(){
    return inputFileName;
}

int Defaults::setInputFileName(string valueIn){
    inputFileName = valueIn;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}

bool Defaults::isStdIoInput(){
    if( inputFileName== "" )
        return true;
    else 
        return false;
}


bool Defaults::getListDevices(){
    return listDevices;
}

int Defaults::setListDevices(bool valueIn){
    listDevices = valueIn;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}

int Defaults::getDevice(){
    return device;
}

int Defaults::setDevice(int valueIn){
    device = valueIn;
    
    return CAMPAIGN_DEFAULTS_SUCCESS;
}


/*
 int main(int argc, const char **argv){
 
 
 Defaults myDefaults(argc, argv);
 
 
 Parameters myParameters;
 
 if( myParameters.tokenizeCommandLine(argc, argv) == PARAMETERS_ERROR){
 std::cout << "did not tokenize correctly! " << std::endl;
 exit(1);
 } else {
 myParameters.display();
 }
 
 std::cout << "testing params" << std::endl;
 
 
 return 0;
 }
 
 */
