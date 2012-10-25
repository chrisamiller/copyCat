#!/bin/bash
version=$(grep "Version" readDepth/DESCRIPTION | awk '{print $2}')
R CMD build readDepth && R CMD INSTALL readDepth_$version.tar.gz
