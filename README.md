# PH_Domain_PIP_Binding_Prediction-
System Requirements: Any operating system where Matlab R2016a or higher version can be installed.
Software requirements: The script was written and tested in Matlab R2016a. The Matlab R2016a was installed on a windows 7 or a windows 10 system. This script should be supported by any operating system that supports Matlab R2016a or higher versions. The script doesn't require any non-standard hardware.
Installation Guide: Please follow standard installation instruction for Matlab installation.
Expected Output: The script will produce two excel files. One for the RFC matrix, another for the SRFC score. The SRFC score Excel file contains three different sheets with SRFC values for SiMPull verified PIP binding proteins, non-binding proteins, and PH domains that have not been experimentally tested.
Instructions: The script requires a text file containing aligned PH domain sequences in fasta fromat. The fasta file contains the PH domain sequences in the following order: 1) PIP binding sequences (1-n), 2) PIP non-binding sequences (n+1 - m), 3) PH domain not experimentally tested (m+1 - end)
