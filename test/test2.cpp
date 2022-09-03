#include <iostream>
#include <omp.h>   // NEW ADD

using namespace std;

int main()
{
    // #pragma omp parallel for num_threads(4) // NEW ADD
    // for(int i=0; i<10; i++)
    // {
    // cout << i << endl;
    // }
    int a = 10;
    cout << a << endl;
    a += 10;
    cout << a << endl;
    a = a + -5;
    cout << a << endl;
    return 0;
}