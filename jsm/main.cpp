#include <stdio.h>
#include <time.h>
#include <random>
#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>

using namespace std;

#define SIGMA 0.100 //0.7071

void singlestep(float *agents, int size, float time_step, int rds){
    random_device rd; //获取随机数种子
    mt19937 gen(rd() + time(0) + rds); //Standard mersenne_twister_engine seeded with rd()
    normal_distribution<double> nd(0,1);
    float bar_v = 0;
    for (int i = 0; i < size; i++){
        bar_v += agents[i];
    }
    bar_v = bar_v/(size - 1);
    // cout << bar_v << endl;

    for (int i = 0; i < size; i++){
        float factor = 1 - size * time_step / (size - 1) + SIGMA * nd(gen) * sqrtf(2 * time_step);
        agents[i] = agents[i] * factor + bar_v * time_step;

        // // // J \prop i
        // float J = 2 * i / (size * (size - 1));
        // float factor = 1 - J * size * time_step / (size - 1) + SIGMA * nd(gen) * sqrtf(2 * time_step);
        // agents[i] = agents[i] * factor + J * bar_v * time_step;

        // // J = -1
        // float factor = 1 + size * time_step / (size - 1) + SIGMA * nd(gen) * sqrtf(2 * time_step);
        // agents[i] = agents[i] * factor - bar_v * time_step;

        // pure noise
        // float factor = 1 + SIGMA * nd(gen) * sqrtf(2 * time_step);
        // agents[i] = agents[i] * factor;
    }

    return ;
}

void singlestep_ex(float *agents, int size, float time_step, int rds){
    random_device rd; //获取随机数种子
    mt19937 gen(rd() + time(0) + rds); //Standard mersenne_twister_engine seeded with rd()
    normal_distribution<double> nd(0,1);
    uniform_int_distribution<int> dist(0,size - 1);
    uniform_real_distribution<float> dist_f(0, 1);

    float J = 0.1;
    float lambda = 0.7;
    for (int i = 0; i < size; i++){
        int a1 = dist(gen);
        int a2 = dist(gen);
        if (a1 == a2) continue;
        // if ( agents[a1]<0 || agents[a2]<0 ) continue;
        // float exch = J * (agents[a1] + agents[a2]);
        // float exch = 0.3 * ( agents[a2]) * dist_f(gen);
        float exch = (lambda - 1) * agents[a1] + (1 - lambda) * ( agents[a1] + agents[a2]) * dist_f(gen);
        // if (agents[a2] < exch) continue;
        agents[a1] += exch;
        agents[a2] -= exch;
    }
    return ;
}

void singlestep_noise(float *agents, int size, float time_step, int rds){
    random_device rd; //获取随机数种子
    mt19937 gen(rd() + time(0) + rds); //Standard mersenne_twister_engine seeded with rd()
    normal_distribution<double> nd(0,1);
    uniform_real_distribution<float> dist(-1,1);

    // float bar_v = 0;
    // for (int i = 0; i < size; i++){
    //     bar_v += agents[i];
    // }
    // bar_v = bar_v/(size - 1);
    // cout << bar_v << endl;

    float J = 0.1;
    float* dv = new float[size];
    for (int i = 0; i < size; i++) dv[i] = 0;

    for (int i = 0; i < size; i++){
        for ( int j = i + 1; j < size; j++){
            if ((agents[j]) < 0 || (agents[i] < 0)) continue;
            float exch = J * (agents[j] + agents[i]) * SIGMA *  (nd(gen) > 0 ? 1 : -1) * sqrtf(2 * time_step);
            // float exch = J * (agents[j] + agents[i]) * SIGMA *  (nd(gen)) * sqrtf(2 * time_step);

            // float exch = J * (agents[j] + agents[i]) * SIGMA * sqrtf(2 * time_step) * (dist(gen) < (agents[i] - agents[j])/(agents[i] + agents[j]) ? 1 : -1);
            dv[i] += exch;
            dv[j] -= exch;
        }
    }

    for (int i = 0; i < size; i++) agents[i] += dv[i];

    delete[] dv;
    return ;
}



float cal_vi_m(float *agents, int size){
    float vim = 0;
    for (int i = 0; i < size; i++){
        vim += agents[i];
    }
    return vim / size;
}

float cal_vi2(float *agents, int size){
    float vi2_m = 0;
    for (int i = 0; i < size; i++){
        vi2_m += agents[i] * agents[i];
    }
    return vi2_m / size;
}

float cal_vivj(float *agents, int size){
    double vivj = 0;
    int num = 0;
    for (int i = 0; i < size; i++){
        for (int j = i + 1; j < size; j++){
            // vivj += agents[i] * agents[j];
            // num ++;
            num ++;
            vivj = vivj + (agents[i] * agents[j] - vivj) / num;
        }
    }
    return vivj;
}

void cal_mean_all_noise(float* vm_noise, float* vims, int T, int rep_t){
    // float vm_noise[T] = {0};
    for (int i = 0; i < rep_t; i++){
        for ( int t = 0; t < T; t++){
            vm_noise[t] += vims[i*T + t] /rep_t ;
        }
    }
}

float cal_mean(float* numbers, int size){
    float ns = 0;
    for ( int i = 0; i < size; i++){
        ns += numbers[i];
    }
    return ns / size;
}

void save_snap(string path, float *agents, int size, float repeat_times){
    ofstream write_f(path);
    for (int t = 0; t < repeat_times; t++){
        for (int i = 0; i < size; i++){
            write_f << agents[i+t*size] << ",";
        }
        write_f << endl;
    }
}

int noiseaver(){
    clock_t start = clock(), end;
    string datapath("./data.txt");
    // string datapath("./data_time30.txt");
    // string datapath("./data_no_corr_mn_11.txt");
    ofstream write_f(datapath);
    random_device rd; //获取随机数种子
    // mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // normal_distribution<double> nd(0,1);
    const int size_agents = 1000;
    const float delta_T = 0.001;
    const int Times = 1000/delta_T;
    const int repeat_times = 1000;
    // const float sigma = 0.7071; // sigma^2 = 0.5
    float* agents_rep = new float[size_agents*repeat_times]; // 数组的每一行是一个 realization
    for ( int i = 0; i < size_agents*repeat_times; i++) agents_rep[i] = 1;

    float vim = 0, vi2 = 0, vivj = 0;
    float* vims = new float[repeat_times];
    float* vi2s = new float[repeat_times];
    float* vivjs = new float[repeat_times];
    for ( int t = 0; t < Times; t++){
        #pragma omp parallel for num_threads(20)
        for ( int realiza = 0; realiza < repeat_times; realiza++){ // 对不同的 realization 计算
            if ( int(realiza%100) == 0) cout << realiza << "\t";
            vims[realiza] = cal_vi_m( agents_rep+size_agents*realiza, size_agents);
            vi2s[realiza] = cal_vi2( agents_rep+size_agents*realiza, size_agents);
            vivjs[realiza] = cal_vivj( agents_rep+size_agents*realiza, size_agents);
            singlestep( agents_rep+size_agents*realiza, size_agents, delta_T, rd());
            // singlestep_noise( agents_rep+size_agents*realiza, size_agents, delta_T, rd());
            // singlestep_ex( agents_rep+size_agents*realiza, size_agents, delta_T, rd());
        }
        vim = cal_mean(vims, repeat_times);
        vi2 = cal_mean(vi2s, repeat_times);
        vivj = cal_mean(vivjs, repeat_times);
        write_f << t << "," << vim << "," << vi2 << "," << vivj << endl;

        if ( int(t%100) == 0) {
            cout << "\n\n" << t << "\n" << endl;
            string a = "./data."+to_string(t)+".txt";
            save_snap(a, agents_rep, size_agents, repeat_times);
        }
    }

    end = clock() - start;
    cout << "time use (s): \t" << (double) end / CLOCKS_PER_SEC;

    delete[] agents_rep;
    delete[] vims;
    delete[] vi2s;
    delete[] vivjs;

    return 0;
}

// 把所有时间都储存在数组中
int noiseall(){
    string datapath("./data.txt");
    ofstream write_f(datapath);
    random_device rd; //获取随机数种子
    // mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // normal_distribution<double> nd(0,1);
    const int size_agents = 1;
    const float delta_T = 0.01;
    const int Times = 1000/delta_T;
    const int repeat_times = 1000;
    // float vims[repeat_times*Times];
    // float vi2s[repeat_times*Times];
    // float vivjs[repeat_times*Times];
    float *vims = new float[repeat_times*Times];
    float *vi2s = new float[repeat_times*Times];
    float *vivjs = new float[repeat_times*Times];
    cout << repeat_times << endl;

    for (int realiza = 0; realiza < repeat_times; realiza++){
        cout << realiza << "\t";
        float agents[size_agents];
        for (int i =0; i < size_agents; i++) agents[i] = 1; //初始化为 1

        float vim, vi2, vivj;

        for (int ts = 0; ts < Times; ts++){
            vim = cal_vi_m(agents, size_agents);
            vi2 = cal_vi2(agents, size_agents);
            vivj = cal_vivj(agents, size_agents);
            // cout << vim << "\t" << vi2 << "\t" << vivj << "\t" << endl;
            vims[Times*realiza + ts] = vim;
            vi2s[Times*realiza + ts] = vi2;
            vivjs[Times*realiza + ts] = vivj;
            // singlestep(agents, size_agents, delta_T, realiza);
            singlestep(agents, size_agents, delta_T, rd());
        }
    }
    cout << "\n" << endl;

    float vmm[Times] = {0};
    cal_mean_all_noise( vmm, vims, Times, repeat_times);
    float vi2m[Times] = {0};
    cal_mean_all_noise( vi2m, vi2s, Times, repeat_times);
    float vivjm[Times] = {0};
    cal_mean_all_noise( vivjm, vivjs, Times, repeat_times);
    for (int ts = 0; ts < Times; ts++){
        // cout << vmm[ts] << "\t";
        write_f << ts << "," << vmm[ts] << "," << vi2m[ts] << "," << vivjm[ts] << endl;
    }
    return 0;
}

int main(){
    noiseaver();
    // noiseall();
}
