#!/bin/bash

echo " Setting Watchman Environment "

# if not set in .bashrc then set here
# export WM_PARENT=/path/to/
 
# testing code
export WM_CODE=${WM_PARENT}Wavedump_Wrapper/

export WM_DATA=${WM_CODE}Data_Storage/
export WM_BINARY=${WM_DATA}Binary_Data/
export WM_ROOT=${WM_DATA}/Raw_Root_Data/

export WM_CONVERT=${WM_CODE}Data_Processing/Binary_Conversion/
export WM_COOK=${WM_CODE}Data_Processing/Cooking/
export WM_COMMON=${WM_CODE}/Common_Tools/
export WM_ANALYSE=${WM_CODE}Data_Analysis/

export CALIBRATION=${WM_ANALYSE}Calibration/
export BIN_TO_ROOT=${WM_ANALYSE}BinToRoot/
export GAIN_TEST_DIR=${WM_ANALYSE}Gain_Test/
export SHIPPING_DATA=${WM_ANALYSE}Shipping_Data/

# header files
export CPATH=${CPATH}:${WM_COMMON}
export CPATH=${CPATH}:${WM_COOK}

# binaries
export PATH=${PATH}:${WM_CONVERT}
export PATH=${PATH}:${WM_CONVERT}/DT_Version/

if [[ "$OSTYPE" == "linux-gnu" ]]; then
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${WM_COOK}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${WM_COMMON}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${WM_ANALYSE}PMT_Analysis
fi

nominal_HV(){
    echo $(grep " $1 " ${WM_COMMON}HVScan.txt) | cut -d " " -f 7
}
