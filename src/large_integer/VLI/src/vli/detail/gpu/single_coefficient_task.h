 #ifndef COEFF_TASK_H
 #define COEFF_TASK_H


namespace vli {
    namespace detail {

	struct single_coefficient_task {
		unsigned short step_count; // how many time the coeff wil be calculate
		unsigned char output_degree_x;
		unsigned char output_degree_y;
		unsigned char output_degree_z;
		unsigned char output_degree_w;
	};

        bool single_coefficient_task_sort(single_coefficient_task const & i, single_coefficient_task const & j);

    } // end namespace detail
} // end name space vli
#endif
