#pragma once

namespace maquis
{
	struct single_coefficient_task
	{
		unsigned short step_count; // how many time the coeff wil be calculate
		unsigned char output_degree_x;
		unsigned char output_degree_y;
	};

	bool single_coefficient_task_sort(single_coefficient_task i, single_coefficient_task j);
}
