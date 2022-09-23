#!/bin/bash

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

inputimage="${1:-examples/dog06.png}"
outputpts="${2:-examples/dog06.pts}"

echo "Image: $inputimage";
echo "Out file: $outputpts"

case "$(uname -s)" in

   Darwin)
     echo 'Mac OS X'
     ;;

   Linux)
     echo 'Linux'
     python prediction/usemodel.py --model prediction/best_model_checkpoint.pth --input $inputimage --output $outputpts
     echo 'Done!'
     ;;

   CYGWIN*|MINGW32*|MSYS*|MINGW*)
     echo 'MS Windows'
     # run python prediction
     python usemodel.py --input $inputimage --output $outputpts
     # execute polyvector with arguments ../my_inputs/alligator_cut2.png ../my_inputs/pythontest.pts
     
     ;;


   *)
     echo 'Other OS' 
     ;;
esac
