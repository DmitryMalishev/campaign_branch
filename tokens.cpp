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
 * \file tokens.cpp
 * \brief Command line tokenization and conversion.
 *
 * \author Author: Marc Sosnick, Contributors: Kai J. Kohlhoff, William Hsu
 * \date 5/26/10
 * \version 1.0
 */

#include "tokens.h"
using namespace std;


Tokens::Tokens(){
}

Tokens::Tokens(int argc, const char **argv){
    tokenizeCommandLine(argc, argv);
}

Tokens::~Tokens(){
}

string Tokens::getArgument(int argNum){
    
    int i;
    
    // iterate through multimap until we reach desired element
    for(paramsIter = params.begin(), i=0; paramsIter != params.end() && i<argNum; ++paramsIter, ++i);
    
    if( paramsIter == params.end() )
        return("");
    else 
        return(paramsIter -> first);
}


int Tokens::getIntegerArgument(string name){
    
    return getIntegerArgument(name, 0);
}


int Tokens::getIntegerArgument(string name, int instance){
    int returnValue;
    string stringToConvert;
    
    if(instance+1 > getArgumentInstances(name)){
        // argument instance is past the number of arguments
        return 0;
    } 
    
    stringToConvert = getStringArgument(name, instance);
    if(stringToConvert == ""){
        // specified argument instance is blank
        return 0;
    }
    
    std::istringstream inpStream(paramsIter->second);
    if( inpStream >> returnValue){
        return returnValue;
    } else {
        // could not convert argument
        return 0;
    }
    
    // should never get here, if we do return error
    return 0;
}

float Tokens::getFloatArgument(string name){
    
    return getFloatArgument(name, 0);
}

float Tokens::getFloatArgument(string name, int instance){
    
    float returnValue;
    string stringToConvert;
    
    if(instance+1 > getArgumentInstances(name)){
        // argument instance is past the number of arguments
        return 0.0;
    } 
    
    stringToConvert = getStringArgument(name, instance);
    if(stringToConvert == ""){
        // specified argument instance is blank
        return 0.0;
    }
    
    std::istringstream inpStream(paramsIter->second);
    if( inpStream >> returnValue){
        return returnValue;
    } else {
        // could not convert argument
        return 0.0;
    }
    
    // should never get here, if we do return error
    return 0.0;
}

string Tokens::getStringArgument(string name){
    return getStringArgument(name, 0);
}


string Tokens::getStringArgument(string name, int instance){
    int i;
    
    if( params.count(name) < instance+1 ){
        return "";
    }
    
    pair<paramsType::iterator, paramsType::iterator> paramsRange;
    paramsType::iterator paramsRangeIter;
    
    paramsIter = params.find(name);
    if( paramsIter != params.end() ){
        paramsRange = params.equal_range(name);
        // iterate until we hit the instance-th parameter
        for(paramsRangeIter = paramsRange.first, i=0; i < instance ; ++paramsRangeIter, ++i);
        return paramsRangeIter->second;
    } else {
        // could not find argument
        return "";
    }
    
    return "";
}

bool Tokens::hasArgument(string paramName){
    
    paramsIter = params.find(paramName);
    if( paramsIter == params.end() ){
        // there is no parameter with that name, 
        // return false.
        return false;
    } else {
        for(int i=0;i < getArgumentInstances(paramName); i++){
            if(paramsIter->second != ""){
                // if the argument is not blank
                return true;
            }
        }
        // all arguments were blank
        return false;
    }
}


bool Tokens::argumentIsNumber(string name){
    
    return argumentIsNumber(name, 0);
}

bool Tokens::argumentIsNumber(string name, int instance){
    int convertValue;
    string stringToConvert;
    
    if( !hasArgument(name) ){
        return false;
    }
    
    if(instance+1 > getArgumentInstances(name)){
        // argument instance is past the number of arguments
        return false;
    } 
    
    stringToConvert = getStringArgument(name, instance);
    if(stringToConvert == ""){
        // specified argument instance is blank
        return false;
    }
    
    // floats will be truncated to ints
    std::istringstream inpStream(paramsIter->second);
    if( inpStream >> convertValue){
        return true;
    } else {
        // could not convert argument
        return false;
    }
    
    // should never get here, if we do return error
    return false;
}

int Tokens::addParameter(string name, string argument){
    params.insert( make_pair(name,argument) );
    return TOKENS_SUCCESS;
}

int Tokens::addParameter(string name){
    params.insert( make_pair(name,""));
    
    return TOKENS_SUCCESS;
}

int Tokens::getArgumentInstances(string name){
    
    return params.count(name);
    
}

int Tokens::getParameterCount(){
    return params.size();
}

void Tokens::display(){
    for(paramsIter = params.begin(); paramsIter != params.end(); ++paramsIter){
        std::cout << "params[" << paramsIter->first << "] = " << paramsIter->second << std::endl;
    }
    
}


int Tokens::tokenizeCommandLine(int argc, const char **argv){
    string s;
    
    // current argument is not an argument, continue processing rest of arguments
    // any arguments unpaired with parameters will be in the _UNPAIRED list
    for(int i=1 ; i<argc; i++){
        s.assign(argv[i]);
        
        if( s.find('-') == 0){
            // if first character is a dash, then it is a parameter
            int dashCount = 0;
            string deDashedArgument = s;
            
            while( deDashedArgument.find('-') == 0 ){
                dashCount++;
                if(dashCount > 2){
                    std::cout << "Error: unknown argument " << s << std::endl;
                    return TOKENS_FAILURE;
                }
                deDashedArgument.assign(deDashedArgument, 1, deDashedArgument.length());
            }
            
            // is there a next argument? if so, see if its a parameter
            if(argc -1 > i) {
                // there is a next argument, test it to see if its an argument
                string paramTestString;
                paramTestString.assign(argv[i+1]);
                if( paramTestString.find('-') == 0 ){
                    // next argv is a parameter, push this parameter
                    addParameter(deDashedArgument);
                } else {
                    // next argv was not a parameter.  push with next arg.
                    addParameter(deDashedArgument, paramTestString);
                    ++i; // increment i to make up for use of the next string
                }
            }  else {
                // no more arguments.  push parameter
                addParameter(deDashedArgument);
            }
        } else {
            // first character is not a dash, add to _UNPAIRED list
            addParameter("_UNPAIRED",s);
        }
    }
    return TOKENS_SUCCESS;
}

