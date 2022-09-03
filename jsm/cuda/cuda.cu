#include <cuda_runtime.h>
#include <stdio.h>
#include <time.h>
#include <random>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

using namespace std;

#define SIGMA 0.7071

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif
#if defined(_WIN32)
#include <chrono>
int gettimeofday(struct timeval* tp, struct timezone* tzp) {
    namespace sc = std::chrono;
    sc::system_clock::duration d = sc::system_clock::now().time_since_epoch();
    sc::seconds s = sc::duration_cast<sc::seconds>(d);
    tp->tv_sec = s.count();
    tp->tv_usec = sc::duration_cast<sc::microseconds>(d - s).count();

    return 0;
}
#endif // _WIN32


double cpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return((double)tp.tv_sec+(double)tp.tv_usec*1e-6);
}

#define CHECK(call)\
{\
    const cudaError_t error=call;\
    if(error!=cudaSuccess)\
    {\
        printf("ERROR: %s:%d,",__FILE__,__LINE__);\
        printf("code:%d,reason:%s\n",error,cudaGetErrorString(error));\
        exit(1);\
    }\
}

int getThreadNum()
{
    cudaDeviceProp prop;
    int count;

    // HANDLE_ERROR(cudaGetDeviceCount(&count));
    cudaGetDeviceCount(&count);
    printf("gpu num %d\n", count);
    // HANDLE_ERROR(cudaGetDeviceProperties(&prop, 0));
    cudaGetDeviceProperties(&prop, 0);
    printf("max thread num: %d\n", prop.maxThreadsPerBlock);
    printf("max grid dimensions: %d, %d, %d)\n",
            prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
    return prop.maxThreadsPerBlock;
}

void initDevice(int devNum)
{
    int dev = devNum;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp,dev));
    printf("Using device %d: %s\n",dev,deviceProp.name);
    CHECK(cudaSetDevice(dev));
}

__global__ void singlerealization(float *agents, int size, int repeat_times, float time_step, float* rds ){
    // int realiza = threadIdx+blockDim*blockIdx;
    int realiza = threadIdx.x+blockDim.x*blockIdx.x;
    if ( realiza < repeat_times ){
        // if (realiza%1024 == 0) printf("!!%d\n",realiza);
        // random_device rd; //获取随机数种子
        // mt19937 gen(rd() + time(0) + rds); //Standard mersenne_twister_engine seeded with rd()
        // normal_distribution<double> nd(1,1);
        float bar_v = 0;
        for (int i = 0; i < size; i++){
            bar_v += agents[i+size*realiza];
        }
        bar_v = bar_v/(size - 1);
        // cout << bar_v << endl;

        for (int i = 0; i < size; i++){
            float factor = 1 - size * time_step / (size - 1) + SIGMA * rds[i+size*realiza] * sqrtf(2 * time_step);
            agents[i+size*realiza] = agents[i+size*realiza] * factor + bar_v * time_step;
            // agents[i+size*realiza] = agents[i+size*realiza] + 1;
            // float factor = 1 + SIGMA * rds[i+size*realiza] * sqrtf(2 * time_step);
            // agents[i+size*realiza] = agents[i+size*realiza] * factor;
            // agents[i+size*realiza] = rds[i+size*realiza];
            // agents[i+size*realiza] = i+size*realiza + rds[i+size*realiza];
        }
    }
}

__global__ void cal_vi_m(float *agents, int size, int repeat_times, float* vims_dev ){
    // int realiza = threadIdx+blockDim*blockIdx;
    int realiza = threadIdx.x+blockDim.x*blockIdx.x;
    if ( realiza < repeat_times ){
        float vim = 0;
        for (int i = 0; i < size; i++){
            vim += agents[i+size*realiza];
        }
        vims_dev[realiza] = vim / size;
    }
}

__global__ void cal_vi2(float *agents, int size, int repeat_times, float* vi2s_dev ){
    // int realiza = threadIdx+blockDim*blockIdx;
    int realiza = threadIdx.x+blockDim.x*blockIdx.x;
    if ( realiza < repeat_times ){
        float vim = 0;
        for (int i = 0; i < size; i++){
            vim += agents[i+size*realiza] * agents[i+size*realiza];
        }
        vi2s_dev[realiza] = vim / size;
    }
}

__global__ void cal_vivj(float *agents, int size, int repeat_times, float* vivjs_dev ){
    // int realiza = threadIdx+blockDim*blockIdx;
    int realiza = threadIdx.x+blockDim.x*blockIdx.x;
    if ( realiza < repeat_times ){
        double vivj = 0;
        int num = 0;
        for (int i = 0; i < size; i++){
            for (int j = i + 1; j < size; j++){
                // vivj += agents[i+size*realiza] * agents[j+size*realiza];
                num ++;
                // vivj = agents[i+size*realiza] * agents[j+size*realiza] / num + vivj/num+(num-1);
                vivj = vivj + (agents[i+size*realiza] * agents[j+size*realiza] - vivj) / num;
            }
        }
        vivjs_dev[realiza] = vivj;
    }
}

float cal_mean(float* numbes, int size){
    float m = 0;
    for (int i = 0; i < size; i++){
        m += numbes[i];
    }
    return m / size;
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
        // float factor = 1 - size * time_step / (size - 1) + SIGMA * nd(gen) * sqrtf(2 * time_step);
        // agents[i] = agents[i] * factor + bar_v * time_step;
        float factor = 0 + SIGMA * nd(gen) * sqrtf(2 * time_step);
        agents[i] = agents[i] * factor;
    }

    return ;
}

int main(){
    initDevice(0);

    random_device rd; //获取随机数种子
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // uniform_real_distribution<float> distr(-1, 1);
    normal_distribution<float> nd(0,1);
    string savepath = "./cuda_data_time001.txt";
    ofstream write_f(savepath);
    // for (int i = 0; i < 5; i++){
    //     cout << nd(gen) << "\t";
    // }
    // getThreadNum();
    const int size_agents = 10000;
    const float delta_T = 0.01;
    const int Times = 1000/delta_T;
    const int repeat_times = 10240;
    cout << repeat_times << endl;

    int nBytes = sizeof(float)*repeat_times;
    float* agents_rep_host = (float*)malloc(nBytes*size_agents);
    float* agents_rep_from_dev = (float*)malloc(nBytes*size_agents);
    float* noise_rds_host = (float*)malloc(nBytes*size_agents);
    float* vims_host = (float*)malloc(nBytes);
    float* vi2s_host = (float*)malloc(nBytes);
    float* vivjs_host = (float*)malloc(nBytes);
    for (int i = 0; i < repeat_times*size_agents; i++) agents_rep_host[i] = 1;

    // ofstream data("noise.txt");
    // for (int i = 0; i < repeat_times*size_agents; i++) data << noise_rds_host[i] << ",";


    // return 0;
    // cudaMalloc, 开辟 device 内存
    float* agents_rep_dev = NULL ;
    float* noise_rds_dev = NULL ;
    float* vims_dev = NULL ;
    float* vi2s_dev = NULL ;
    float* vivjs_dev = NULL ;

    CHECK(cudaMalloc((void**)&agents_rep_dev, nBytes*size_agents)); 
    CHECK(cudaMalloc((void**)&noise_rds_dev, nBytes*size_agents)); 
    CHECK(cudaMalloc((void**)&vims_dev, nBytes));
    CHECK(cudaMalloc((void**)&vi2s_dev, nBytes));
    CHECK(cudaMalloc((void**)&vivjs_dev, nBytes));

    //输入数据从主机内存拷贝到设备内存
    CHECK(cudaMemcpy(agents_rep_dev, agents_rep_host, nBytes*size_agents, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(noise_rds_dev, noise_rds_host, nBytes*size_agents, cudaMemcpyHostToDevice));
    // CHECK(cudaMemcpy(vims_dev, vims_host, nBytes, cudaMemcpyHostToDevice));
    // CHECK(cudaMemcpy(vi2s_dev, vi2s_host, nBytes, cudaMemcpyHostToDevice));
    // CHECK(cudaMemcpy(vivjs_dev, vivjs_host, nBytes, cudaMemcpyHostToDevice));

    string a = "./cuda_data_time001.-1.txt";
    save_snap(a, agents_rep_host, size_agents, repeat_times);

    // //一维线程块，32×32
    dim3 block(1024);
    // //一维线程网格
    dim3 grid((repeat_times-1)/block.x+1);
    cout << block.x << "\t" << grid.x << endl;
    double gpuStart = cpuSecond();
    cout << "start gpu" << endl;
    for (int t = 0; t < Times; t++) {
        for (int i = 0; i < repeat_times*size_agents; i++){
            noise_rds_host[i] = nd(gen);
        }
        CHECK(cudaMemcpy(noise_rds_dev, noise_rds_host, nBytes*size_agents, cudaMemcpyHostToDevice));

        singlerealization<<<grid, block >>>( agents_rep_dev, size_agents, repeat_times, delta_T, noise_rds_dev);
        cal_vi_m<<<grid, block >>>( agents_rep_dev, size_agents, repeat_times, vims_dev);
        cal_vi2<<<grid, block >>>( agents_rep_dev, size_agents, repeat_times, vi2s_dev);
        cal_vivj<<<grid, block >>>( agents_rep_dev, size_agents, repeat_times, vivjs_dev);
        CHECK(cudaDeviceSynchronize());

        CHECK(cudaMemcpy(vims_host, vims_dev, nBytes, cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(vi2s_host, vi2s_dev, nBytes, cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(vivjs_host, vivjs_dev, nBytes, cudaMemcpyDeviceToHost));
        float vim1;
        float vim2;
        float vim3;
        vim1 = cal_mean(vims_host, repeat_times);
        vim2 = cal_mean(vi2s_host, repeat_times);
        vim3 = cal_mean(vivjs_host, repeat_times);
        if (t%100 == 0) cout << t << "\t";
        // write_f << t << "," << vim1 << "," << endl;
        write_f << t << "," << vim1 << "," << vim2 << "," << vim3 << endl;

        if (t%1000 == 0){
            CHECK(cudaMemcpy(agents_rep_host, agents_rep_dev, nBytes*size_agents, cudaMemcpyDeviceToHost));
            string a = "./cuda_data_time001."+to_string(t)+".txt";
            save_snap(a, agents_rep_host, size_agents, repeat_times);
        }
    }

    double gpuTime = cpuSecond() - gpuStart;
    printf("GPU Execution Time: %f sec\n", (double) gpuTime);


    double cpuStart=cpuSecond();
    for ( int realiza = 0; realiza < repeat_times; realiza++){ // 对不同的 realization 计算
        // if ( int(realiza%1000) == 0) cout << realiza << "\t";
        singlestep( agents_rep_host+size_agents*realiza, size_agents, delta_T, rd());
    }
    double cpuTime = cpuSecond() - cpuStart;
    printf("CPU Execution Time: %f sec\n", (double) cpuTime);

    CHECK(cudaMemcpy(agents_rep_from_dev, agents_rep_dev, nBytes*size_agents, cudaMemcpyDeviceToHost));
    // CHECK(cudaMemcpy(vims_host, vims_dev, nBytes, cudaMemcpyDeviceToHost));
    // CHECK(cudaMemcpy(vi2s_host, vi2s_dev, nBytes, cudaMemcpyDeviceToHost));
    // CHECK(cudaMemcpy(vivjs_host, vivjs_dev, nBytes, cudaMemcpyDeviceToHost));

    // for (int i = 0; i < 10; i++){
    //     cout << agents_rep_from_dev[i] << "\t" << agents_rep_host[i] << "\t";
    // }

    cudaFree(agents_rep_dev);
    free(agents_rep_from_dev);
    free(agents_rep_host);
    cudaDeviceReset();

    return 0;
}