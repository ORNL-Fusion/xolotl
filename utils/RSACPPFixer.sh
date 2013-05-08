#!/bin/bash

#------------------------------------------------------------------------------
#
# File: RSACPPFixer.sh
# Author: Jay Jay Billings
# Author Contact: billingsjj@ornl.gov
# Description: This script uses sed, indent and astyle to fixed deemed problems 
#              in C++ files generate from IBM RSA 8.x UML models. This script
#              should only be used on files that are *freshly* generated from
#              from IBM RSA.
# Date: 20130410
#
#------------------------------------------------------------------------------

# Kill @generated tags
find . -type f -print0 | xargs -0 sed -i '\/\/@generated "UML to C++ (com.ibm.xtools.transform.uml2.cpp.CPPTransformation)"/d'
# Kill <p>
find . -type f -print0 | xargs -0 sed -i 's/<p>/ /' 
# Kill </p>
find . -type f -print0 | xargs -0 sed -i 's/<\/p>/ /'
# Kill "Begin" blocks. It looks for lines containing //Begin section, //TODO: 
# Add or End section and removes them. 
find . -type f -print0 | xargs -0 sed -i '/\/\/Begin section/d;/\/\/TODO: Add/d;/\/\/End section/d'
# Kill multiple empty lines and replace them with only one.
find . -type f -print0 | xargs -0 sed -i '/^$/N;/^\n$/D'
