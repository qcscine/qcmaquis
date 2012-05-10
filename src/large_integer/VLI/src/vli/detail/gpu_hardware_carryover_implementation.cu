#include <vector>
#include <algorithm>

#include "gpu_hardware_carryover_implementation.h"
#include "single_coefficient_task.h"
#include "compile_time_constants.h"
#include "numeric.h"
#include "kernels_gpu_asm.hpp"



__device__ vli::detail::single_coefficient_task execution_plan[MUL_BLOCK_SIZE * MAX_ITERATION_COUNT];
__device__ unsigned int workblock_count_by_warp[MUL_BLOCK_SIZE / 32];

// (low_32_bits : acc, high_32_bits : co) = acc + c1 * c2
#define __wide_mad(acc,co,c1,c2) asm("{mad.lo.cc.u32 %0, %2, %3, %0;\n	madc.hi.u32 %1, %2, %3, 0;}" : "+r"(acc), "=r"(co) : "r"(c1), "r"(c2));

// (low_32_bits : acc, high_32_bits : co) = acc + c1 * c2 + co
#define __wide_mad_with_carry_in(acc,co,c1,c2) asm("{\n	.reg .u32 nv;\n	mad.lo.cc.u32 nv, %2, %3, %1;\n	madc.hi.u32 %1, %2, %3, 0;\n	add.cc.u32 %0, nv, %0;\n	addc.u32 %1, %1, 0;}" : "+r"(acc), "+r"(co) : "r"(c1), "r"(c2));

__global__ void
__launch_bounds__(MUL_BLOCK_SIZE, 2)
polynomial_mul_full(
	const unsigned int * __restrict__ in1,
	const unsigned int * __restrict__ in2,
	const unsigned int element_count,
	unsigned int * __restrict__ out)
{
	__shared__ unsigned int in_buffer1[SLICE_PADDED * DEGREE_BOUND_Y * INT_DEGREE];
	__shared__ unsigned int in_buffer2[SLICE_PADDED * DEGREE_BOUND_Y * INT_DEGREE];

	const unsigned int local_thread_id = threadIdx.x;
        const unsigned int element_id = blockIdx.x;

	// Copy both input polynomials into the shared memory
	{
		const unsigned int * in_shifted1 = in1 + (element_id * (DEGREE_BOUND_X * DEGREE_BOUND_Y * INT_DEGREE));
		const unsigned int * in_shifted2 = in2 + (element_id * (DEGREE_BOUND_X * DEGREE_BOUND_Y * INT_DEGREE));
		unsigned int index_id = local_thread_id;
		#pragma unroll
		for(unsigned int i = 0; i < (DEGREE_BOUND_X * DEGREE_BOUND_Y * INT_DEGREE) / MUL_BLOCK_SIZE; ++i)
		{
			unsigned int coefficient_id = index_id / INT_DEGREE;
			unsigned int degree_id = index_id % INT_DEGREE;
			unsigned int current_degree_y = coefficient_id / DEGREE_BOUND_X;
			unsigned int current_degree_x = coefficient_id % DEGREE_BOUND_X;
			unsigned int local_index_id = current_degree_x + (current_degree_y * SLICE_PADDED) + (degree_id * (SLICE_PADDED * DEGREE_BOUND_Y));
			in_buffer1[local_index_id] = in_shifted1[index_id];
			in_buffer2[local_index_id] = in_shifted2[index_id];
			index_id += MUL_BLOCK_SIZE;
		}
		if (index_id < (DEGREE_BOUND_X * DEGREE_BOUND_Y * INT_DEGREE))
		{
			unsigned int coefficient_id = index_id / INT_DEGREE;
			unsigned int degree_id = index_id % INT_DEGREE;
			unsigned int current_degree_y = coefficient_id / DEGREE_BOUND_X;
			unsigned int current_degree_x = coefficient_id % DEGREE_BOUND_X;
			unsigned int local_index_id = current_degree_x + (current_degree_y * SLICE_PADDED) + (degree_id * (SLICE_PADDED * DEGREE_BOUND_Y));
			in_buffer1[local_index_id] = in_shifted1[index_id];
			in_buffer2[local_index_id] = in_shifted2[index_id];
		}
		__syncthreads();
	}

        unsigned int c1[INT_DEGREE],c2[INT_DEGREE];
	unsigned int res[2*INT_DEGREE];
	unsigned int res1[2*INT_DEGREE];

	unsigned int iteration_count = workblock_count_by_warp[local_thread_id / 32];
	for(unsigned int iteration_id = 0; iteration_id < iteration_count; ++iteration_id)
	{
		vli::detail::single_coefficient_task task = execution_plan[local_thread_id + (iteration_id * MUL_BLOCK_SIZE)];
		const unsigned int step_count = task.step_count;

		if (step_count > 0)
		{
			const unsigned int output_degree_y = task.output_degree_y;
			const unsigned int output_degree_x = task.output_degree_x;

			#pragma unroll
			for(unsigned int i = 0; i < 2*INT_DEGREE; ++i)
				res[i] = 0;

			const unsigned int start_degree_x_inclusive = output_degree_x > (DEGREE_BOUND_X - 1) ? output_degree_x - (DEGREE_BOUND_X - 1) : 0;
			const unsigned int end_degree_x_inclusive = output_degree_x < DEGREE_BOUND_X ? output_degree_x : (DEGREE_BOUND_X - 1);
			unsigned int current_degree_x = start_degree_x_inclusive;
			unsigned int current_degree_y = output_degree_y > (DEGREE_BOUND_Y - 1) ? output_degree_y - (DEGREE_BOUND_Y - 1) : 0;
			for(unsigned int step_id = 0; step_id < step_count; ++step_id) {
                            unsigned int * in_polynomial1 = in_buffer1 + current_degree_x + (current_degree_y * SLICE_PADDED);
                            unsigned int * in_polynomial2 = in_buffer2 + (output_degree_x - current_degree_x) + ((output_degree_y - current_degree_y) * SLICE_PADDED);

                            #pragma unroll
                            for(unsigned int i = 0; i < 2*INT_DEGREE; ++i)
                               res1[i] = 0;

                            #pragma unroll
                            for(unsigned int degree1 = 0; degree1 < INT_DEGREE; ++degree1)
                                c1[degree1] = in_polynomial1[degree1 * (SLICE_PADDED * DEGREE_BOUND_Y)];
 
                            #pragma unroll
                            for(unsigned int degree2 = 0; degree2 < INT_DEGREE; ++degree2)
                                c2[degree2] = in_polynomial2[degree2  * (SLICE_PADDED * DEGREE_BOUND_Y)];
 
                            unsigned int sign = (c1[INT_DEGREE-1]>>31) ^ (c2[INT_DEGREE-1]>>31);

                            if(c1[INT_DEGREE-1] >> 31 != 0)
                                vli::detail::negate192_gpu(c1); 

                            if(c2[INT_DEGREE-1] >> 31 != 0)
                                vli::detail::negate192_gpu(c2); 

                            vli::detail::mul384_384_gpu(res1,c1,c2);

                            if(sign != 0)
                               vli::detail::negate384_gpu(res1);

                            vli::detail::add384_384_gpu(res,res1);
                                   
				// Calculate the next pair of input coefficients to be multiplied and added to the result
                             current_degree_x++;
                             if (current_degree_x > end_degree_x_inclusive) {
                                 current_degree_x = start_degree_x_inclusive;
                                 current_degree_y++;
                             }
			}

			unsigned int coefficient_id = output_degree_y * MULT_RESULT_DEGREE_BOUND_X + output_degree_x;
			unsigned int * out2 = out + (coefficient_id * element_count *2* INT_DEGREE) + element_id; // coefficient->int_degree->element_id
			#pragma unroll
			for(unsigned int i = 0; i < 2*INT_DEGREE; ++i) {
				// This is a strongly compute-bound kernel,
				// so it is fine to waste memory bandwidth by using non-coalesced writes in order to have less instructions,
				//     less synchronization points, less shared memory used (and thus greater occupancy) and greater scalability.
				*out2 = res[i];
				out2 += element_count;
			}
		} // if (step_count > 0)
	} //for(unsigned int iteration_id
}

__global__ void
__launch_bounds__(SUM_BLOCK_SIZE, 2)
polynomial_sum_intermediate_full(
	const unsigned int * __restrict__ intermediate,
	const unsigned int element_count,
	unsigned int * __restrict__ out)
{
	__shared__ unsigned int buf[SUM_BLOCK_SIZE * 2*INT_DEGREE_PADDED];

	unsigned int local_thread_id = threadIdx.x;
	unsigned int coefficient_id = blockIdx.x;

	unsigned int * t1 = buf + (local_thread_id * 2*INT_DEGREE_PADDED);
	#pragma unroll
	for(unsigned int i = 0; i <2* INT_DEGREE; ++i)
		t1[i] = 0;

	const unsigned int * in2 = intermediate + (coefficient_id * element_count *2* INT_DEGREE) + local_thread_id;
	for(unsigned int element_id = local_thread_id; element_id < element_count; element_id += SUM_BLOCK_SIZE)
	{
		unsigned int carry_over = 0;
		#pragma unroll
		for(unsigned int degree = 0; degree < 2*INT_DEGREE; ++degree)
		{
			unsigned long long res_wide = (unsigned long long)t1[degree] + in2[degree * element_count] + carry_over;
			t1[degree] = res_wide;
			carry_over = res_wide >> 32;
		}
		in2 += SUM_BLOCK_SIZE;
	}

	#pragma unroll
	for(unsigned int stride = SUM_BLOCK_SIZE >> 1; stride > 0; stride >>= 1)
	{
		__syncthreads();

		if (local_thread_id < stride)
		{
			unsigned int * t2 = buf + ((local_thread_id + stride) *2* INT_DEGREE_PADDED);

			unsigned int carry_over = 0;
			#pragma unroll
			for(unsigned int degree = 0; degree < 2*INT_DEGREE; ++degree)
			{
				unsigned long long res_wide = (unsigned long long)t1[degree] + t2[degree] + carry_over;
				t1[degree] = res_wide;
				carry_over = res_wide >> 32;
			}
		}
	}

	if (local_thread_id == 0)
	{
		unsigned int * out2 = out + (coefficient_id * 2*INT_DEGREE);
		#pragma unroll
		for(unsigned int i = 0; i < 2*INT_DEGREE; ++i)
			out2[i] = buf[i];
	}
}

namespace vli {
    namespace detail {
	gpu_hardware_carryover_implementation::gpu_hardware_carryover_implementation()
		: d_intermediate_result(0), element_count_prepared(0)
	{
		std::vector<unsigned int> workblock_count_by_warp_local(MUL_BLOCK_SIZE / 32);
		for(unsigned int i = 0; i < workblock_count_by_warp_local.size(); ++i)
			workblock_count_by_warp_local[i] = 0;
		std::vector<unsigned int> work_total_by_size(MUL_BLOCK_SIZE / 32);
		for(unsigned int i = 0; i < work_total_by_size.size(); ++i)
			work_total_by_size[i] = 0;

		std::vector<vli::detail::single_coefficient_task> tasks(((MULT_RESULT_DEGREE_BOUND_X*MULT_RESULT_DEGREE_BOUND_Y + 32 - 1) / 32) * 32);
		for(unsigned int degree_y = 0; degree_y < MULT_RESULT_DEGREE_BOUND_Y; ++degree_y)
		{
			for(unsigned int degree_x = 0; degree_x < MULT_RESULT_DEGREE_BOUND_X; ++degree_x)
			{
				vli::detail::single_coefficient_task& task = tasks[degree_y * MULT_RESULT_DEGREE_BOUND_X + degree_x];
				task.output_degree_x = degree_x;
				task.output_degree_y = degree_y;

				task.step_count = (std::min<unsigned int>((MULT_RESULT_DEGREE_BOUND_X - 1) - degree_x, degree_x) + 1) * (std::min<unsigned int>((MULT_RESULT_DEGREE_BOUND_Y - 1) - degree_y, degree_y) + 1);
			}
		}
		// Fill the task list up to the multiple of the warp size
		for(unsigned int i = MULT_RESULT_DEGREE_BOUND_X * MULT_RESULT_DEGREE_BOUND_Y; i < tasks.size(); ++i)
		{
			vli::detail::single_coefficient_task& task = tasks[i];
			task.output_degree_x = 0;
			task.output_degree_y = 0;
			task.step_count = 0;
		}
		// Sort the tasks in step_count descending order
		std::sort(tasks.begin(), tasks.end(), single_coefficient_task_sort);

		vli::detail::single_coefficient_task empty_task;
		empty_task.step_count = 0;
		std::vector<vli::detail::single_coefficient_task> tasks_reordered(MUL_BLOCK_SIZE * MAX_ITERATION_COUNT, empty_task);
		for(unsigned int batch_id = 0; batch_id < tasks.size() / 32; ++batch_id)
		{
			unsigned int warp_id = std::min_element(work_total_by_size.begin(), work_total_by_size.end()) - work_total_by_size.begin();
			std::copy(
				tasks.begin() + (batch_id * 32),
				tasks.begin() + ((batch_id + 1) * 32),
				tasks_reordered.begin() + (workblock_count_by_warp_local[warp_id] * MUL_BLOCK_SIZE) + (warp_id * 32));
			unsigned int max_step_count = tasks[batch_id * 32].step_count;

			workblock_count_by_warp_local[warp_id]++;
			work_total_by_size[warp_id] += max_step_count;
		}

		cudaMemcpyToSymbolAsync(workblock_count_by_warp, &(*workblock_count_by_warp_local.begin()), sizeof(unsigned int) * workblock_count_by_warp_local.size());
		cudaMemcpyToSymbolAsync(execution_plan, &(*tasks_reordered.begin()), sizeof(vli::detail::single_coefficient_task) * tasks_reordered.size());
	}

	gpu_hardware_carryover_implementation::~gpu_hardware_carryover_implementation()
	{
//		if (d_intermediate_result)
//			cudaFree(d_intermediate_result);
	}

	void  gpu_hardware_carryover_implementation::run(
		unsigned int * d_input_vector1,
		unsigned int * d_input_vector2,
		unsigned int * d_output_polynomial,
		unsigned int element_count)
	{
		prepare(element_count);

		{
			dim3 grid(element_count);
			dim3 threads(MUL_BLOCK_SIZE);

			polynomial_mul_full<<<grid,threads>>>(
				d_input_vector1,
				d_input_vector2,
				element_count,
				d_intermediate_result);

//			getLastCudaError("CUDA polynomial_mul Kernel execution failed");
		}

		{
			dim3 grid(MULT_RESULT_DEGREE_BOUND_X * MULT_RESULT_DEGREE_BOUND_Y);
			dim3 threads(SUM_BLOCK_SIZE);

			polynomial_sum_intermediate_full<<<grid,threads>>>(
				d_intermediate_result,
				element_count,
				d_output_polynomial);

//			getLastCudaError("CUDA polynomial_sum_intermediate Kernel execution failed");
		}

//		return cudaSuccess;
	}

	void gpu_hardware_carryover_implementation::prepare(unsigned int max_element_count)
	{
		if (max_element_count > element_count_prepared)
		{
			if (d_intermediate_result)
			    (cudaFree(d_intermediate_result));
			(cudaMalloc((void **)&d_intermediate_result, 2*SINGLE_MULT_RESULT_POLYNOMIAL_UINT_COUNT * sizeof(unsigned int) * max_element_count));
		}
	}
    }
}
