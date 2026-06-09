to run it you need to do it:


g++ -O2 -o femto main.cpp $(root-config --cflags --libs) -lMathCore -lGenVector


./femto input.txt
