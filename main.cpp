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

using namespace std;

int main(int argc, char *argv[])
{
    int x = 10;
    int *y = new int[5];
    cout << "hello world" << endl;
    for (int i = 0; i < argc; i++)
    {
        int *yy = new int[5];
        int zz = i * 20;
        cout << argv[i] << endl;
    }
    return 0;
}
