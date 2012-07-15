#include "single_coefficient_task.h"


namespace vli {
    namespace detail {

        bool single_coefficient_task_sort(single_coefficient_task const & i, single_coefficient_task const &j){
		if (i.step_count != j.step_count)
			return (i.step_count > j.step_count);

		if (i.output_degree_w != j.output_degree_w)
			return (i.output_degree_w < j.output_degree_w);

		if (i.output_degree_z != j.output_degree_z)
			return (i.output_degree_z < j.output_degree_z);

		if (i.output_degree_y != j.output_degree_y)
			return (i.output_degree_y < j.output_degree_y);

		return (i.output_degree_x < j.output_degree_x);
	}
    } // end namespace detail
} // end name space vli

