all:
	clang++ -std=c++20 -fopenmp -Ofast *.cc -lomp
