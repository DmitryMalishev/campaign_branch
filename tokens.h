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
 * \file tokens.h
 * \brief Header file for command line tokenization and conversion.
 *
 * \author Author: Marc Sosnick, Contributors: Kai J. Kohlhoff, William Hsu
 * \date 5/26/10
 * \version 1.0
 */
#ifndef __TOKENS_H__
#define __TOKENS_H__

#include <string>
#include <map>
#include <iostream>
#include <sstream>

//using namespace std;

/** @def TOKENS_SUCCESS 
 *  @brief Return value upon success
 **/
#define TOKENS_SUCCESS  0

/** @def TOKENS_FAILURE
 *  @brief Return value upon error
 **/
#define TOKENS_FAILURE  1

/**
 * @brief Tokenizes command line and provides argument testing.
 *
 * Tokenizes command line into parameter-argument pairs.  Allows unlimited
 * multiple pairings of arguments to parameters.  Provides type testing
 * and conversion for arguments.  All arguments are stored as the original
 * std::string internally, until the specific type is requested for that argument.<br>
 * Arguments can be passed unpaired, e.g. multiple filenames.  Any such unpaired
 * arguments are paired with the _UNPAIRED parameter name, and can thus be
 * retreived in order.<br>
 * Parameters may have a maximum of two dashes (i.e. - or --, not --- or ----).
 * Leading dashes are *not* stored in params, e.g. -N is equivalent to --N. 
 * 
 * @author Marc Sosnick
 * @date 5/26/10
 **/
class Tokens{
    typedef std::multimap<std::string, std::string> paramsType;
    paramsType params;
    paramsType::iterator paramsIter;
    
public:
    Tokens();
    ~Tokens();
    
    /**
     * @brief Tokenizes passed command line parameters.
     **/
    Tokens(int argc, const char **argv);
    
    /**
     * @brief Gets the specified parameter's specified argument as an integer.
     * Gets the specified parameter's specified argument as an integer.  Use 
     * in conjunction with hasArgument() and argumentIsNumber(), and 
     * getArgumentInstances() to determine if parameterName does have argument 
     * if it is a number, and how many arguments.
     * @param name - name of the parameter from which to get the argument.
     * @param instance - if more than one argument exists for that parameter, 
     *                       returns the insatnce-th argument. 
     * @return the specified parameter's argument as an integer
     * If there is no argument associated with name, returns 0.  If there is no 
     * instance-th argument, returns 0.  
     **/
    int getIntegerArgument(std::string name, int instance);
    
    /**
     * @brief Gets the specified parameter's first argument as an integer.
     * @param name - name of the parameter to get the argument from
     * @return the specified parameter's argument as an integer
     **/
    int getIntegerArgument(std::string name);
    
    /**
     * @brief Gets the specified parameter's specified argument as a float.
     * Gets the specified parameter's specified argument as a float.  Use 
     * in conjunction with hasArgument() and argumentIsNumber(), and 
     * getArgumentInstances() to determine if parameterName does have argument 
     * if it is a number, and how many arguments.
     * @param name - name of the parameter from which to get the argument.
     * @param instance - if more than one argument exists for that parameter, 
     *                       returns the insatnce-th argument. 
     * @return the specified parameter's argument as a float. 
     * If there is no argument associated with name, returns 0.0.  If there is no 
     * instance-th argument, returns 0.0.  
     **/
    float   getFloatArgument(std::string name, int instance);
    
    /**
     * @brief Gets the specified parameter's first argument as an float.
     * @param name - name of the parameter to get the argument from
     * @return the specified parameter's argument as a float 
     **/
    float   getFloatArgument(std::string name);
    
    /**
     * @brief Gets the specified parameter's specified argument as a std::string.
     * @param name - name of the parameter from which to get argument.
     * @param instance - if more than one argument exists for that parameter, 
     *                       returns the insatnce-th argument. 
     * @return the specified parameter's argument as a std::string. 
     * If there is no argument associated with name, returns empty std::string.  If there is no 
     * instance-th argument, returns empty std::string.  
     **/
    std::string  getStringArgument(std::string name, int instance);
    
    /**
     * @brief Gets the specified parameter's first argument as an std::string.
     * @param name - name of the parameter to get the argument from
     * @return the specified parameter's argument as a float 
     **/
    std::string  getStringArgument(std::string name);
    
    /**
     * @brief Indicates if the specified parameter has a non-blank argument.
     * @return If name has an argument that is not blank, returns true, 
     * otherwise returns false. If name is not a parameter, returns false.
     **/
    bool    hasArgument(std::string parameter);
    
    /**
     * @brief Returns the number of parameters, including duplicates.
     * @return the number of parameters stored, including non-unique
     **/
    int     getParameterCount();
    
    /**
     * @brief Returns the number of arguments associated with specified parameter.
     * @return number of arguments
     **/
    int getArgumentInstances(std::string name);
    
    /**
     * @brief Returns the argumentNumber-th argument as a std::string.
     * @param argument number.  
     * @return Returns argumentNumber-th argument as a std::string.  If argumentNumber is 
     *         larger than the number of elements, returns empty std::string.
     **/
    std::string  getArgument(int argumentNumber);
    
    /** 
     * @brief Indicates if specified parameter's first argument is a number.
     * @return If name's fist argument is a number, returns true, otherwise
     * returns false.  If name is not a parameter, returns false.
     **/
    bool    argumentIsNumber(std::string name);
    
    /**
     * @brief Indicates if specified parameter's specified argument is a number.
     * @return If name's instance-th argument is a number, returns true, otherwise
     * returns false.  If name is not a parameter, returns false.
     **/
    bool    argumentIsNumber(std::string name, int instance);
    
    /**
     * @brief Adds specified parameter and argument to the list of parameters.
     * @return TOKENS_FAILURE if there is an error, otherwise returns TOKENS_SUCCESS.
     **/
    int addParameter(std::string name, std::string argument);
    
    /**
     * @brief Adds specified parameter without argument to the list of parameters.
     * @return TOKENS_FAILURE if there is an error, otherwise returns TOKENS_SUCCESS.
     **/
    int addParameter(std::string name);
    
    /**
     * @brief Tokenizes the command line.
     * @return TOKENS_FAILURE if there is an error, otherwise returns TOKENS_SUCCESS.
     **/
    int     tokenizeCommandLine(int argc, const char **argv);
    
    /**
     * @brief Displays all tokenized parameters and their values.
     **/    
    void    display();
};

#endif



