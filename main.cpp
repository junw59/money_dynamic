#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <set>
#include <omp.h>
#include <chrono>

using namespace std;

int gettimeofday(struct timeval* tp, struct timezone* tzp) {
    namespace sc = std::chrono;
    sc::system_clock::duration d = sc::system_clock::now().time_since_epoch();
    sc::seconds s = sc::duration_cast<sc::seconds>(d);
    tp->tv_sec = s.count();
    tp->tv_usec = sc::duration_cast<sc::microseconds>(d - s).count();

    return 0;
}

double cpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return((double)tp.tv_sec+(double)tp.tv_usec*1e-6);
}

int main(int argc, char *argv[])
{
    double duration;
    // clock_t start, finish;
    // start = clock();
    double start = cpuSecond();
    int x = 10;
    int *y = new int[5];
    cout << "hello world" << endl;
    for (int i = 0; i < argc; i++)
    {
        int *yy = new int[5];
        int zz = i * 20;
        cout << i << argv[i] << endl;
    }


    double last_time = cpuSecond() - start;
    cout << last_time << endl;
    // finish = clock();
    // duration = (double)(finish - start);
    // cout << duration/CLOCKS_PER_SEC << endl;
    return 0;
}
