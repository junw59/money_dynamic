#include <stdio.h>
#include <time.h>
#include <random>
#include <iostream>
#include <fstream>
#include <windows.h>

using namespace std;

void genn(float a, int j){
    random_device rd; //获取随机数种子
	mt19937 gen(j); //Standard mersenne_twister_engine seeded with rd()
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
	mt19937 gen(rd()+time(NULL)); //Standard mersenne_twister_engine seeded with rd()
	// uniform_real_distribution<float> dist(-1, 1);
    uniform_int_distribution<int> dist(0,400 - 1);
    uniform_real_distribution<float> dist_f(0, 1);
    // normal_distribution<float> nd(0,1);
    int size = 10000;
    // float* agents = new float[size*10];

    // for (int i = 0; i < size*10; i++) agents[i] = i;

    // for (int i = 0; i < 10; i++){
    //     out_array(agents + i * 10, size);
    //     // cout << agents[i] << "\t";
    // }
    // Sleep(500);
    // end = clock() - start;
    // cout << "times \t" << (float) end/CLOCKS_PER_SEC;
    // ofstream fileout("data_real.txt");
    // float mea = 0;
    for (int i = 0; i < 500; i++){
        // fileout << dist(gen) << endl;
        // mea += nd(gen);
        cout << dist_f(gen) << "\t";
        // cout << rd() << "\t";
        // genn(rd(), i);
    }
    // cout << mea / 1000000;

    // float* arra = new float[20];
    // for ( int i = 0; i < 20; i++){
    //     arra[i] = i;
    // }
    // for ( int i = 0; i < 20; i++){
    //     cout << (arra + i)[0] << "\t";
    // }
    // float vivj = 0;
    // int num = 0;
    // for (int i = 0; i < size; i++){
    //     for (int j = i + 1; j < size; j++){
    //         // vivj += 1;
    //         num ++;
    //         // vivj = 1 / num + vivj/num*(num-1);
    //         vivj = (1 - vivj) / num + vivj;
    //         // vivj += 1;
    //     }
    // }

    // cout << vivj << "\t" << vivj / num << "\t" << num << endl;
    return 0;
}
