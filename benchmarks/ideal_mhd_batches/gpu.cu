#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <vector>
#include "ideal_mhd_batches.cuh"
extern "C" void mhd_cpu_reference(int family, int n, const double* x, double* y);
static void checked(cudaError_t rc) {
    if (rc != cudaSuccess) throw std::runtime_error(cudaGetErrorString(rc));
}
template<int family, int mode, typename Input, typename Output>
static void run(int n, int repeats) {
    constexpr int counts[4]={20,27,7,23};
    constexpr int m=counts[family-1], ni=2*m+1, no=m+2;
    std::vector<double> x(std::size_t(n)*ni,0), reference(std::size_t(n)*no);
    auto put=[&](int i,int column,double value){x[std::size_t(column-1)*n+i]=value;};
    for (int i=0;i<n;++i) {
        const double u=0.25+((i+1)%97)/128.0, v=-0.125+((i+1)%53)/64.0;
        put(i,ni,1);
        if constexpr (family <= 2) {
            put(i,1,1); put(i,7,i%2 ? -1 : 1); put(i,8,1); put(i,9,1);
            if constexpr (family == 1) {put(i,18,u); put(i,20,v);}
            else {put(i,19,u);put(i,21,v);put(i,27,2);put(i,15,3);put(i,22,0.5);}
        }
        if constexpr (family == 3) {
            put(i,1,i%2 ? -1 : 1);put(i,2,2);put(i,3,4);put(i,5,3);
            put(i,6,u);put(i,7,v);
        }
        if constexpr (family == 4) {put(i,3,1);put(i,5,1);put(i,12,u);put(i,21,v);}
        for (int k=1;k<=m;++k) put(i,m+k,(k%2 ? -1.0 : 1.0)/64);
    }
    mhd_cpu_reference(family,n,x.data(),reference.data());
    std::vector<Input> input(x.begin(),x.end());
    std::vector<Output> output(std::size_t(n)*no+16,Output(-987));
    Input* dx; Output* dy;
    checked(cudaMalloc(&dx,input.size()*sizeof(Input)));
    checked(cudaMalloc(&dy,output.size()*sizeof(Output)));
    cudaEvent_t begin,end;
    checked(cudaEventCreate(&begin)); checked(cudaEventCreate(&end));
    auto launch=[&]() {
        if constexpr (family == 1 && mode == 0)
            mhd_two_fp64_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 1 && mode == 1)
            mhd_two_fp32_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 1 && mode == 2)
            mhd_two_mixed_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 2 && mode == 0)
            mhd_three_fp64_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 2 && mode == 1)
            mhd_three_fp32_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 2 && mode == 2)
            mhd_three_mixed_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 3 && mode == 0)
            mhd_mass_fp64_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 3 && mode == 1)
            mhd_mass_fp32_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 3 && mode == 2)
            mhd_mass_mixed_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 4 && mode == 0)
            mhd_xi_fp64_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 4 && mode == 1)
            mhd_xi_fp32_batch<<<(n+127)/128,128>>>(n,dx,dy);
        if constexpr (family == 4 && mode == 2)
            mhd_xi_mixed_batch<<<(n+127)/128,128>>>(n,dx,dy);
        checked(cudaGetLastError());
    };
    double max_error=0, scale=0, checksum=0;
    auto validate=[&]() {
        checked(cudaMemcpy(output.data(),dy,output.size()*sizeof(Output),cudaMemcpyDeviceToHost));
        max_error=0;scale=0;checksum=0;
        for (std::size_t i=0;i<reference.size();++i) {
            if (!std::isfinite(output[i])) throw std::runtime_error("nonfinite CUDA MHD output");
            max_error=std::max(max_error,std::abs(double(output[i])-reference[i]));
            scale=std::max(scale,std::abs(reference[i]));
        }
        const double budget=mode == 0 ? 2e-12 : 512*std::numeric_limits<float>::epsilon();
        if (max_error > budget*(1+scale)) throw std::runtime_error("CUDA MHD precision oracle");
        for (std::size_t i=reference.size();i<output.size();++i)
            if (output[i] != Output(-987)) throw std::runtime_error("CUDA MHD tail overwrite");
        for (int i=0;i<n;++i) checksum+=double(output[i]); // FP64 accumulation
    };
    for (int sample=-1;sample<=repeats;++sample) {
        // Transfer-inclusive and resident times are separately reported.
        checked(cudaEventRecord(begin));
        checked(cudaMemcpyAsync(dx,input.data(),input.size()*sizeof(Input),cudaMemcpyHostToDevice));
        checked(cudaMemcpyAsync(dy,output.data(),output.size()*sizeof(Output),cudaMemcpyHostToDevice));
        launch();
        checked(cudaMemcpyAsync(output.data(),dy,output.size()*sizeof(Output),cudaMemcpyDeviceToHost));
        checked(cudaEventRecord(end));checked(cudaEventSynchronize(end));
        float transfer_ms;
        checked(cudaEventElapsedTime(&transfer_ms,begin,end));
        checked(cudaEventRecord(begin));
        launch();
        checked(cudaEventRecord(end));checked(cudaEventSynchronize(end));
        float resident_ms;
        checked(cudaEventElapsedTime(&resident_ms,begin,end));
        validate();
        if (sample>0) std::printf("%d,%d,%d,%d,%.9g,%.9g,%.9g,%.17g,%zu\n",
            family,mode,n,sample,resident_ms*1e-3,transfer_ms*1e-3,max_error,checksum,
            input.size()*sizeof(Input)+output.size()*sizeof(Output));
    }
    checked(cudaEventDestroy(begin));checked(cudaEventDestroy(end));
    checked(cudaFree(dx));checked(cudaFree(dy));
}
int main(int argc,char** argv) {
    const int n=argc>1 ? std::atoi(argv[1]) : 257, repeats=argc>2 ? std::atoi(argv[2]) : 3;
    if (n<1 || repeats<1) return 2;
    cudaDeviceProp device; checked(cudaGetDeviceProperties(&device,0));
    std::fprintf(stderr,"GPU: %s; arch %d%d\n",device.name,device.major,device.minor);
    run<1,0,double,double>(n,repeats);
    run<1,1,float,float>(n,repeats);
    run<1,2,float,double>(n,repeats);
    run<2,0,double,double>(n,repeats);
    run<2,1,float,float>(n,repeats);
    run<2,2,float,double>(n,repeats);
    run<3,0,double,double>(n,repeats);
    run<3,1,float,float>(n,repeats);
    run<3,2,float,double>(n,repeats);
    run<4,0,double,double>(n,repeats);
    run<4,1,float,float>(n,repeats);
    run<4,2,float,double>(n,repeats);
}
