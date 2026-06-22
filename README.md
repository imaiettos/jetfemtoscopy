to run it you need to do it:


g++ -O2 -o femto main.cpp $(root-config --cflags --libs) -lMathCore -lGenVector


./femto input.txt input_wFSI.txt

input.txt contains "test_data.root"


input_wFSI.txt contains "pp_parton_cascade_batch0_1.root"
