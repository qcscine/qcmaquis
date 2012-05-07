#pragma once

#include <cuda_runtime.h>

namespace vli {
    namespace detail {

	class gpu_hardware_carryover_implementation {
	public:
		gpu_hardware_carryover_implementation();

		virtual ~gpu_hardware_carryover_implementation();

		void prepare(unsigned int max_element_count);

		// All buffers should be in device memory
		void run( unsigned int * d_input_vector1, unsigned int * d_input_vector2, unsigned int * d_output_polynomial, unsigned int element_count);

	private:
		// Disable copying
		gpu_hardware_carryover_implementation(const gpu_hardware_carryover_implementation&);
		gpu_hardware_carryover_implementation& operator =(const gpu_hardware_carryover_implementation&);

		unsigned int element_count_prepared;
		unsigned int * d_intermediate_result;
	};
    } // end namespace detail
 }//end namespace vli

