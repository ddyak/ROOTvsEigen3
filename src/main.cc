


#include <benchmark/benchmark.h>

// Define benchmark
static void BM_SomeFoo(benchmark::State& state) {
    for (auto _ : state) {
        foo(); // performance measured code
    }
}
// Register the function as a benchmark
BENCHMARK(BM_StringCreation);

BENCHMARK_MAIN();

//g++ main.cpp -std=c++11 -lbenchmark -lpthread -O2 -o benchmark
// 3, 4, 6, 7

#include <benchmark/benchmark.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <KineTypedefs.h>
#include <KineTools.h>
#include <TLorentzVector.h>

#define BENCH "notmatterword"
#include <iostream>
#include <vector>

#include <Eigen/Dense>

Matrix3f m;
m << 1, 2, 3,
     4, 5, 6,
     7, 8, 9;
std::cout << m;

Output:
1 2 3
4 5 6
7 8 9

std::vector<double> arr;
#ifdef BENCH

static void BM_ThreeVectorCreation(benchmark::State &state)
{
  for (auto _ : state)
  {
    sct::ThreeVector a(1, 2, 3);
    benchmark::DoNotOptimize(sct::kine::costh(a));
  }
}
BENCHMARK(BM_ThreeVectorCreation);

static void BM_TVector3Creation(benchmark::State &state)
{
  for (auto _ : state)
  {
    TVector3 a(1, 2, 3);
    benchmark::DoNotOptimize(a.CosTheta());
  }
}
BENCHMARK(BM_TVector3Creation);


static void BM_EigenMatrixMult3x3(benchmark::State &state)
{
  const int size = 3;
  for (int i = 0; i < 1000; ++i)
  {
    arr.push_back(rand() - 0.5);
  }
  sct::SquareMatrix<size> A, B;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
        B(i, j) = arr[(++k + 500) % 1000];
      }
    }
    benchmark::DoNotOptimize(A*B);
  }
}
BENCHMARK(BM_EigenMatrixMult3x3);

static void BM_TMatrixDMult3x3(benchmark::State &state)
{
  const int size = 3;
  TMatrixD A(size, size), B(size, size);
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
        B(i, j) = arr[(++k + 500) % 1000];
      }
    }
    A*B;
  }
}
BENCHMARK(BM_TMatrixDMult3x3);

static void BM_EigenMatrixMult4x4(benchmark::State &state)
{
  const int size = 4;
  sct::SquareMatrix<size> A, B;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
        B(i, j) = arr[(++k + 500) % 1000];
      }
    }
    benchmark::DoNotOptimize(A*B);
  }
}
BENCHMARK(BM_EigenMatrixMult4x4);

static void BM_TMatrixDMult4x4(benchmark::State &state)
{
  const int size = 4;
  TMatrixD A(size, size), B(size, size);
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
        B(i, j) = arr[(++k + 500) % 1000];
      }
    }
    benchmark::DoNotOptimize(A*B);
  }
}
BENCHMARK(BM_TMatrixDMult4x4);

static void BM_EigenMatrixMult7x7(benchmark::State &state)
{
  const int size = 7;
  sct::SquareMatrix<size> A, B;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
        B(i, j) = arr[(++k + 500) % 1000];
      }
    }
    benchmark::DoNotOptimize(A*B);
  }
}
BENCHMARK(BM_EigenMatrixMult7x7);

static void BM_TMatrixDMult7x7(benchmark::State &state)
{
  const int size = 7;
  TMatrixD A(size, size), B(size, size);
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
        B(i, j) = arr[(++k + 500) % 1000];
      }
    }
    benchmark::DoNotOptimize(A*B);
  }
}
BENCHMARK(BM_TMatrixDMult7x7);


static void BM_EigenMatrixInv3x3(benchmark::State &state)
{
  const int size = 3;
  sct::SquareMatrix<size> A;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
      }
    }
    benchmark::DoNotOptimize(A.inverse());
  }
}
BENCHMARK(BM_EigenMatrixInv3x3);

static void BM_TMatrixDInv3x3(benchmark::State &state)
{
  const int size = 3;
  TMatrixD A(size, size);
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
      }
    }
    benchmark::DoNotOptimize(A.Invert());
  }
}
BENCHMARK(BM_TMatrixDInv3x3);

static void BM_EigenMatrixInv4x4(benchmark::State &state)
{
  const int size = 4;
  sct::SquareMatrix<size> A;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
      }
    }
    benchmark::DoNotOptimize(A.inverse());
  }
}
BENCHMARK(BM_EigenMatrixInv4x4);

static void BM_TMatrixDInv4x4(benchmark::State &state)
{
  const int size = 4;
  TMatrixD A(size, size);
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
      }
    }
    benchmark::DoNotOptimize(A.Invert());
  }
}
BENCHMARK(BM_TMatrixDInv4x4);


static void BM_EigenMatrixInv7x7(benchmark::State &state)
{
  const int size = 7;
  sct::SquareMatrix<size> A;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
      }
    }
    benchmark::DoNotOptimize(A.inverse());
  }
}
BENCHMARK(BM_EigenMatrixInv7x7);

static void BM_TMatrixDInv7x7(benchmark::State &state)
{
  const int size = 7;
  TMatrixD A(size, size);
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
      }
    }
    benchmark::DoNotOptimize(A.Invert());
  }
}
BENCHMARK(BM_TMatrixDInv7x7);


static void BM_Eigen_Matrix_4x4_vec4(benchmark::State &state)
{
  const int size = 4;
  sct::SquareMatrix<size> A;
  Eigen::Vector4d D;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      D(i) = arr[(k+++500)%1000];
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
      }
    }
    benchmark::DoNotOptimize(A*D);
  }
}
BENCHMARK(BM_Eigen_Matrix_4x4_vec4);

static void BM_Root_Matrix4x4_vec4(benchmark::State &state)
{
  const int size = 4;
  TMatrixD A(size, size), B(size, 1);
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      B(i, 0) = arr[(++k+500) % 1000];
      for (int j = 0; j < size; ++j)
      {
        A(i, j) = arr[++k % 1000];
      }
    }
    benchmark::DoNotOptimize(A*B);
  }
}
BENCHMARK(BM_Root_Matrix4x4_vec4);

static void BM_Eigen_Norm(benchmark::State &state)
{
  const int size = 4;
  sct::SquareMatrix<size> A;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
        A(i) = arr[++k % 1000];
    }
    benchmark::DoNotOptimize(A.norm());
  }
}
BENCHMARK(BM_Eigen_Norm);

static void BM_Root_Norm(benchmark::State &state)
{
  const int size = 4;
  TLorentzVector A;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i) 
    {
        A(i) = arr[++k % 1000];
    }
    benchmark::DoNotOptimize(A.Mag());
  }
}
BENCHMARK(BM_Root_Norm);

static void BM_Eigen_Vec4xVec4(benchmark::State &state)
{
  const int size = 4;
  Eigen::Vector4d A, B;
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      A(i) = arr[(k++ + 500)%1000];
      B(i) = arr[(k) % 1000];
    }
    benchmark::DoNotOptimize(A*B.transpose());
  }
}
BENCHMARK(BM_Eigen_Vec4xVec4);

static void BM_Root_Vec4xVec4(benchmark::State &state)
{
  const int size = 4;
  TMatrixD A(size, 1), B(1, size);
  int k = 0;
  for (auto _ : state)
  {
    for (int i = 0; i < size; ++i)
    {
      B(0, i) = arr[(++k+500) % 1000];
      A(i, 0) = arr[(++k) % 1000];
    }
    benchmark::DoNotOptimize(B*A);
  }
}
BENCHMARK(BM_Root_Vec4xVec4);


BENCHMARK_MAIN();

#else

int main()
{
  for (int i = 0; i < 1000; ++i)
  {
    arr.push_back(rand() - 0.5);
  }
  TMatrixD A(3, 3), B(3, 3);
  for (int k = 0; k < 10000; ++k)
  {
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        A(i, j) = arr[k % 1000];
        B(i, j) = arr[(k + 500) % 1000];
      }
    }

    TMatrixD C(3, 3);
    C.Mult(A, B);
  }
  // C.Print();
  return 0;
}

#endif
