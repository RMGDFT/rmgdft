#include <Tiled-MM/tiled_mm.hpp>
#include <Tiled-MM/device_vector.hpp>
#include <Tiled-MM/util.hpp>
#include <Tiled-MM/gpu_blas_handle.hpp>
#include <Tiled-MM/gpu_blas_api.hpp>
#include <Tiled-MM/gpu_runtime_api.hpp>

#include <iostream>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <random>

#include <cxxopts.hpp>

void compute_reference(double* a, double* b, double* c,
        int m, int n, int k,
	int ld_a, int ld_b, int ld_c,
        double alpha, double beta,
	char trans_a, char trans_b) {
    // sumatrix size to multiply
    int a_m = trans_a == 'N' ? m : k;
    int a_n = trans_a == 'N' ? k : m;

    int b_m = trans_b == 'N' ? k : n;
    int b_n = trans_b == 'N' ? n : k;

    int c_m = m;
    int c_n = n;

    gpu::device_vector<double> a_device(a_m * a_n);
    gpu::device_vector<double> b_device(b_m * b_n);
    gpu::device_vector<double> c_device(c_m * c_n);

    gpu::copy_to_device(a, a_device.data(), a_device.size());
    gpu::copy_to_device(b, b_device.data(), b_device.size());
    gpu::copy_to_device(c, c_device.data(), c_device.size());

    gpu::gpu_blas_handle handle;

    auto transa = gpu::get_blas_operation(trans_a);
    auto transb = gpu::get_blas_operation(trans_b);

    gpu::blas_api::dgemm(handle.handle(), transa, transb,
                         m, n, k, 
			 &alpha, 
			 a_device.data(), ld_a,
                         b_device.data(), ld_b, 
			 &beta, 
			 c_device.data(), ld_c);

    gpu::copy_to_host(c_device.data(), c, c_device.size());
}

using value_type = double;
using size_type  = size_t;

template <typename T>
void fill_matrix(T* ptr, size_t size) {
    static std::random_device dev;                        // seed
    static std::mt19937 rng(dev());                       // generator
    static std::uniform_real_distribution<T> dist(10.0); // distribution

    for (unsigned i = 0; i < size; ++i) {
        ptr[i] = T{dist(rng)};
    }
}

template <typename T>
void copy_matrix(T* from, T* to, std::size_t size) {
    for (unsigned i = 0; i < size; ++i) {
        to[i] = from[i];
    }
}

template <typename T>
bool equal(T* v1, T* v2, size_t len, double eps=1e-6) {
    for (unsigned i = 0; i < len; ++i) {
        auto value1 = *(v1 + i);
        auto value2 = *(v2 + i);
        if (std::abs(value1 - value2) > eps) {
            return false;
        }
    }
    return true;
}

template <typename T>
void print_matrix(T* mat, int m, int n, char trans) {
    if (trans != 'N') {
        std::swap(m, n);
    }

    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            auto el = j * m + i;
            std::cout << mat[el] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "\n";
}

int main(int argc, char** argv){
    cxxopts::Options options("Testing Tiled-MM", "Tests the result of correctness of the Tiled-MM algorithm.");

    options.add_options()
	    ("m,m_dim", 
	         "The number of rows of the resulting matrix C.", 
		 cxxopts::value<int>()->default_value("1000"))
	    ("n,n_dim", 
	         "The number of columns of the resulting matrix C.", 
		 cxxopts::value<int>()->default_value("1000"))
	    ("k,k_dim", 
	         "The size of the shared dimension between matrices A and B.", 
		 cxxopts::value<int>()->default_value("1000"))
	    ("tile_m", 
	        "The tile size for dimension m.", 
		cxxopts::value<int>()->default_value("5000"))
	    ("tile_n", 
	        "The tile size for dimension n.", 
		cxxopts::value<int>()->default_value("5000"))
	    ("tile_k", 
	        "The tile size for dimension k.", 
		cxxopts::value<int>()->default_value("5000"))
	    ("n_streams", 
	        "The number of GPU streams to use.", 
		cxxopts::value<int>()->default_value("2"))
	    ("ld_a", 
	         "The leading dimension of matrix A.", 
		 cxxopts::value<int>()->default_value("0"))
	    ("ld_b", 
	         "The leading dimension of matrix B.", 
		 cxxopts::value<int>()->default_value("0"))
	    ("ld_c", 
	         "The leading dimension of matrix C.", 
		 cxxopts::value<int>()->default_value("0"))
	    ("t,transpose", 
	         "Whether matrices A and B are not-transposed (N), transposed (T) or transpose-conjugated (C). \
	          For example, passing NT means that matrix A is not transposed, whereas the matrix B is transposed.", 
	          cxxopts::value<std::string>()->default_value("NN"))
	    ("alpha", 
	         "The constant alpha in: C = beta*C + alpha*A*B.", 
		 cxxopts::value<double>()->default_value("1.0"))
	    ("beta", "The constant beta in: C = beta*C + alpha*A*B.", 
	         cxxopts::value<double>()->default_value("0.0"))
	    ;

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    auto m = result["m_dim"].as<int>();
    auto n = result["n_dim"].as<int>();
    auto k = result["k_dim"].as<int>();

    auto tile_m = result["tile_m"].as<int>();
    auto tile_n = result["tile_n"].as<int>();
    auto tile_k = result["tile_k"].as<int>();

    auto n_streams = result["n_streams"].as<int>();

    auto transpose = result["transpose"].as<std::string>();
    // transform to upper-case
    std::transform(transpose.begin(), transpose.end(), transpose.begin(),
        [&](char c) {
            return std::toupper(c);
        }
    );
    std::unordered_set<std::string> transpose_options = {
        "NN", "TT", "NT", "TN"
    };

    // check if transpose takes a correct value
    if (std::find(transpose_options.begin(), transpose_options.end(), transpose) == transpose_options.end()) {
        std::cout << "[ERROR]: --transpose option \
        can only take the following values: " << std::endl;
        for (const auto& el : transpose_options) {
            std::cout << el << ", ";
        }
        std::cout << std::endl;
        return 0;
    }

    auto trans_a = transpose[0];
    auto trans_b = transpose[1];

    // matrix sizes to multiply
    int a_m = trans_a == 'N' ? m : k;
    int a_n = trans_a == 'N' ? k : m;

    int b_m = trans_b == 'N' ? k : n;
    int b_n = trans_b == 'N' ? n : k;

    int c_m = m;
    int c_n = n;

    int ld_a = std::max(a_m, result["ld_a"].as<int>());
    int ld_b = std::max(b_m, result["ld_b"].as<int>());
    int ld_c = std::max(c_m, result["ld_c"].as<int>());

    auto alpha = result["alpha"].as<double>();
    auto beta = result["beta"].as<double>();

    int num_reps = 1; // since testing

    bool small_sizes = std::max(m, std::max(n, k)) < 20;

    std::cout << "==================================================" << std::endl;
    std::cout << "                Benchmarking Tiled-MM    " << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << "         MATRIX SIZES " << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << " A = (" << a_m << ", " << a_n << ")" << std::endl;
    std::cout << " B = (" << b_m << ", " << b_n << ")" << std::endl;
    std::cout << " C = (" << c_m << ", " << c_n << ")" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << "         LEADING DIMS " << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << " LD_A = " << ld_a << std::endl;
    std::cout << " LD_B = " << ld_b << std::endl;
    std::cout << " LD_C = " << ld_c << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << "      SCALING CONSTANTS " << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << " alpha = " << alpha << std::endl;
    std::cout << " beta  = " << alpha << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << "      TRANSPOSE FLAGS " << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << " trans_a = " << trans_a << std::endl;
    std::cout << " trans_b = " << trans_b << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << "         TILE SIZES " << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << " tile_m = " << tile_m << std::endl;
    std::cout << " tile_n = " << tile_n << std::endl;
    std::cout << " tile_k = " << tile_k << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << "      ADDITIONAL OPTIONS " << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << " num. of gpu streams = " << n_streams << std::endl;
    std::cout << " num. of repetitions = " << num_reps << std::endl;
    std::cout << "=============================" << std::endl;

    // A dimensions: MxK
    auto a_host = gpu::malloc_pinned<value_type>(a_m * a_n, 1);
    // B dimensions: KxN
    auto b_host = gpu::malloc_pinned<value_type>(b_m * b_n, 1);
    // C dimensions: MxN
    auto c_host = gpu::malloc_pinned<value_type>(c_m * c_n, 0);
    auto c_host2 = gpu::malloc_pinned<value_type>(c_m * c_n, 0);
    auto c_host_reference = gpu::malloc_pinned<value_type>(c_m * c_n, 0);

    fill_matrix(a_host, a_m * a_n);
    fill_matrix(b_host, b_m * b_n);
    fill_matrix(c_host, c_m * c_n);

    copy_matrix(c_host, c_host_reference, c_m * c_n);
    copy_matrix(c_host, c_host2, c_m * c_n);

    if (small_sizes) {
        std::cout << "Initial values in matrix A: " << std::endl;
        print_matrix(a_host, a_m, a_n, trans_a);
        std::cout << "Initial values in matrix B: " << std::endl;
        print_matrix(b_host, b_m, b_n, trans_b);
        std::cout << "Initial values in matrix C: " << std::endl;
        print_matrix(c_host, c_m, c_n, 'N');
    }

    compute_reference(a_host, b_host, c_host_reference, 
		      m, n, k, 
		      ld_a, ld_b, ld_c,
		      alpha, beta, 
		      trans_a, trans_b);

    if (small_sizes) {
        std::cout << "Correct result C = beta*C + alpha*A*B: " << std::endl;
        print_matrix(c_host_reference, m, n, 'N');
    }

    auto ctx = gpu::make_context<double>(n_streams, tile_m, tile_n, tile_k);

    // VERSION WITH COPYING C BACK
    bool copy_c_back = true;
    // compute c = alpha * a * b + beta * c

    auto start = std::chrono::steady_clock::now();
    gpu::gemm(*ctx,
              trans_a, trans_b, 
              m, n, k, 
	      alpha,
	      a_host, ld_a,
	      b_host, ld_b,
	      beta,
	      c_host, ld_c,
	      false, copy_c_back);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Time [ms] with copying C back: " << duration << std::endl;

    if (small_sizes) {
        std::cout << "Computed result by Tiled-MM with copying C back : " << std::endl;
        print_matrix(c_host, m, n, 'N');
    }

    bool correct = equal(c_host, c_host_reference, c_m*c_n);

    // VERSION WITHOUT COPYING C BACK
    // compute the same but don't copy c back
    start = std::chrono::steady_clock::now();
    gpu::gemm(*ctx,
              trans_a, trans_b, 
	      m, n, k, 
	      alpha, 
	      a_host, ld_a,
	      b_host, ld_b,
	      beta, 
	      c_host2, ld_c,
	      false, !copy_c_back);
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Time [ms] without copying C back: " << duration << std::endl;

    gpu::copy_to_host(ctx->get_full_device_buffer_c().data(), c_host2, c_m * c_n);

    if (small_sizes) {
        std::cout << "Computed result by Tiled-MM without copying C back : " << std::endl;
        print_matrix(c_host2, c_m, c_n, 'N');
    }

    correct = correct && equal(c_host2, c_host_reference, c_m*c_n);

    std::cout << "The result is " << (correct ? "CORRECT" : "NOT CORRECT") << std::endl;;

    return !correct;
}


