rm -rf results/

#mpirun -np 6 ../../src/buildVP/hemepure -in input.xml -out results
mpirun -np 5 ../../src/build_WKtest/hemepure -in input.xml -out results

./hemeXtract -X results/Extracted/inlet.dat > results/inlet.txt
./hemeXtract -X results/Extracted/outlet.dat > results/outlet.txt
./hemeXtract -X results/Extracted/planeY.dat > results/planeY.txt

bash toParaview.sh

python3 postProcessFlowRates.py
