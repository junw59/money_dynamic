#include <stdio.h>
#include <time.h>
#include <random>
#include <iostream>
#include <fstream>
#include <windows.h>

using namespace std;

void genn(float a, int j){
    random_device rd; //获取随机数种子
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	// uniform_real_distribution<float> distr(-ETA, ETA);
    normal_distribution<double> nd(0,1);
    for (int i = 0; i < 10; i++){
        // fileout << nd(gen) << endl;
        cout << nd(gen) << "\t";
        // cout << time(0) << "\t";
        // cout << rd() << "\t";
    }
    cout << endl;
}

void out_array(float* agents, int size){
    for ( int i = 0; i < size; i++){
        cout << agents[i] << "\t";
    }
    cout << endl;
}

int main(){
    clock_t start = clock(), end;
    random_device rd; //获取随机数种子
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	// uniform_real_distribution<float> distr(-ETA, ETA);
    normal_distribution<double> nd(0,1);
    int size = 10;
    float* agents = new float[size*10];

    for (int i = 0; i < 5; i++){
        cout << rd() << "\t";
    }

    return 0;
    for (int i = 0; i < size*10; i++) agents[i] = i;

    // for (int i = 0; i < 10; i++){
    //     out_array(agents + i * 10, size);
    //     // cout << agents[i] << "\t";
    // }
    // Sleep(500);
    // end = clock() - start;
    // cout << "times \t" << (float) end/CLOCKS_PER_SEC;
    // ofstream fileout("data.txt");
    // float mea = 0;
    for (int i = 0; i < 5; i++){
        // fileout << nd(gen) << endl;
        // mea += nd(gen);
        // cout << nd(gen) << endl;
        // cout << rd() << "\t";
        genn(rd(), i);
    }
    // cout << mea / 1000000;

    // float* arra = new float[20];
    // for ( int i = 0; i < 20; i++){
    //     arra[i] = i;
    // }
    // for ( int i = 0; i < 20; i++){
    //     cout << (arra + i)[0] << "\t";
    // }
    return 0;
}