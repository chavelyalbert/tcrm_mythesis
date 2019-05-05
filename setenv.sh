#/bin/bash
source activate tcrm201
codePath=`pwd` 
echo $codePath
export PYTHONPATH=$PYTHONPATH:$codePath:$codePath/Utilities
echo $PYTHONPATH
#cd ~/Documents/tcrm-2.0.2
