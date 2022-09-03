#include <iostream>

using namespace std;

extern "C" {
    void prip(int a, int b){
        cout << "hello world!" << a << "\t" << b << endl;
    }

    int add(int a, int b){
        return a + b;
    }

    float cal_vivj(float *agents, int size){
        double vivj = 0;
        int num = 0;
        cout << size << endl;
        for (int i = 0; i < size; i++){
            for (int j = i + 1; j < size; j++){
                // vivj += agents[i] * agents[j];
                // num ++;
                num ++;
                vivj = vivj + (agents[i] * agents[j] - vivj) / num;
                // cout << num << "\t" << vivj << endl;
            }
        }
        return vivj;
    }

    float mean_cij(float *agents, int size, int lins){
        float mv = 0;
        int num = 0;
        for (int i = 0; i < lins; i++){
            num++;
            float vij = cal_vivj(agents + i * size, size);
            mv = mv + (vij - mv) / num;
        }
        return mv;
    }
}
