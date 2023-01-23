#include <Tiled-MM/tiled_mm.hpp>

#include <cxxopts.hpp>

#include <iostream>
#include <cmath>
#include <cstdio>
#include <chrono>

int main(int argc, char** argv){
    cxxopts::Options options("Benchmarking Tiled-MM", "Measures the runtime of the Tiled-MM algorithm.");

    options.add_options()
	    ("m,m_dim", 
	        "the number of rows of the resulting matrix c.", 
		cxxopts::value<int>()->default_value("1000"))
	    ("n,n_dim", 
	        "the number of columns of the resulting matrix c.", 
		cxxopts::value<int>()->default_value("1000"))
	    ("k,k_dim", 
	        "the size of the shared dimension between matrices a and b.", 
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
	    ("r,n_rep", 
	        "The number of repetitions.", 
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
	    ("beta", 
	        "The constant beta in: C = beta*C + alpha*A*B.", 
		cxxopts::value<double>()->default_value("0.0"))
	    ;

    auto result = options.parse(argc, argv);

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
    
    auto num_runs = result["n_rep"].as<int>(); 

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
    std::cout << " num. of repetitions = " << num_runs << std::endl;
    std::cout << "=============================" << std::endl;

    size_t flops_per_mul = 2 * m * n * k;

    // A dimensions: MxK
    auto a_host = gpu::malloc_pinned<double>(a_m * a_n, 1);
    // B dimensions: KxN
    auto b_host = gpu::malloc_pinned<double>(b_m * b_n, 1);
    // C dimensions: MxN
    auto c_host = gpu::malloc_pinned<double>(c_m * c_n, 0);

    auto ctx = gpu::make_context<double>(n_streams, tile_m, tile_n, tile_k);


    std::cout << std::endl;

    std::cout << "==================================================" << std::endl;
    std::cout << "         Results of benchmarking Tiled-MM    " << std::endl;
    std::cout << "==================================================" << std::endl;
    std::cout << " 1) The version with copying C to back to host: " << std::endl;

    auto start = std::chrono::steady_clock::now();
    for(int i = 0; i < num_runs + 1; ++i) {
        // compute c = alpha * a * b + beta * c
        // no need to pin the memory, since matrices already pinned
        if (i == 1) {
            start = std::chrono::steady_clock::now();
        }
        bool pin_buffers = false;
        bool copy_c_back = true;
        gpu::gemm(*ctx, 
		  trans_a, trans_b,
		  m, n, k,
		  alpha,
		  a_host, ld_a,
		  b_host, ld_b,
		  beta,
		  c_host, ld_c,
		  pin_buffers, copy_c_back);
    }

    auto end = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    auto time_per_mul =  time / num_runs;
    std::cout << "    -> Avg Time [ms] = " << time_per_mul << std::endl;

    auto flops = flops_per_mul / (1e-3 * time_per_mul) / 1e9;
    std::cout << "    -> Throughput [Gflops] = " << flops << std::endl;

    std::cout << "==================================================" << std::endl;
    std::cout << " 2) The version without copying C to back to host: " << std::endl;

    start = std::chrono::steady_clock::now();
    for(int i=0; i<num_runs+1; ++i) {
        // compute c = alpha * a * b + beta * c
        // no need to pin the memory, since matrices already pinned
        if (i == 1) {
            start = std::chrono::steady_clock::now();
        }
        bool pin_buffers = false;
        bool copy_c_back = false;
        gpu::gemm(*ctx, 
		  trans_a, trans_b,
		  m, n, k,
		  alpha,
		  a_host, ld_a,
		  b_host, ld_b,
		  beta,
		  c_host, ld_c,
		  pin_buffers, copy_c_back);
    }

    end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    time_per_mul =  time / num_runs;
    std::cout << "    -> Avg Time [ms] = " << time_per_mul << std::endl;

    flops = flops_per_mul / (1e-3 * time_per_mul) / 1e9;
    std::cout << "    -> Throughput [Gflops] = " << flops << std::endl;
    std::cout << "==================================================" << std::endl;

    return 0;
}


