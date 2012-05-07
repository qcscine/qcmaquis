#include "single_coefficient_task.h"

namespace maquis
{
	bool single_coefficient_task_sort(single_coefficient_task i, single_coefficient_task j)
	{
		if (i.step_count != j.step_count)
			return (i.step_count > j.step_count);

		if (i.output_degree_y != j.output_degree_y)
			return (i.output_degree_y < j.output_degree_y);

		return (i.output_degree_x < j.output_degree_x);
	}
}
