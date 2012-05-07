#pragma once

namespace vli {
    namespace detail {

	struct single_coefficient_task {
		unsigned short step_count; // how many time the coeff wil be calculate
		unsigned char output_degree_x;
		unsigned char output_degree_y;
	};

	bool single_coefficient_task_sort(single_coefficient_task i, single_coefficient_task j) {
		if (i.step_count != j.step_count)
			return (i.step_count > j.step_count);

		if (i.output_degree_y != j.output_degree_y)
			return (i.output_degree_y < j.output_degree_y);

		return (i.output_degree_x < j.output_degree_x);
	}
    } // end namespace detail
} // end name space vli
