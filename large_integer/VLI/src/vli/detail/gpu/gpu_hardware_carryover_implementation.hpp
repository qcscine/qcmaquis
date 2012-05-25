
namespace vli {
    namespace detail {

    template <typename BaseInt, std::size_t Size, unsigned int Order>
    gpu_hardware_carryover_implementation<BaseInt, Size, Order>::gpu_hardware_carryover_implementation(){

        unsigned int mul_block_size = (((Order*2-1)*(Order*2-1)/2U >= 256U) ? 256U : (((Order*2-1)*(Order*2-1)/2U+32U-1U)/32U*32U));// 32U is the warp size here
        unsigned int max_iteration_count = ((Order*2-1)*(Order*2-1)+31U)/32U;
        //it's a big constructor
        element_count_prepared=0;
    
        std::vector<unsigned int> workblock_count_by_warp_local(mul_block_size / 32,0);
        std::vector<unsigned int> work_total_by_size(mul_block_size / 32,0);
        std::vector<vli::detail::single_coefficient_task> tasks((((Order*2-1)*(Order*2-1) + 32 - 1) / 32) * 32);
      
        for(unsigned int degree_y = 0; degree_y < (Order*2-1); ++degree_y) {
            for(unsigned int degree_x = 0; degree_x < (Order*2-1); ++degree_x) {
                vli::detail::single_coefficient_task& task = tasks[degree_y * (Order*2-1) + degree_x];
                task.output_degree_x = degree_x;
                task.output_degree_y = degree_y;
                task.step_count = (std::min<unsigned int>(((Order*2-1) - 1) - degree_x, degree_x) + 1) * (std::min<unsigned int>(((Order*2-1) - 1) - degree_y, degree_y) + 1);
            }
        }
        // Fill the task list up to the multiple of the warp size
        for(unsigned int i = (Order*2-1) * (Order*2-1); i < tasks.size(); ++i) {
               vli::detail::single_coefficient_task& task = tasks[i];
               task.output_degree_x = 0;
               task.output_degree_y = 0;
               task.step_count = 0;
        }
       // Sort the tasks in step_count descending order
         std::sort(tasks.begin(), tasks.end(), single_coefficient_task_sort);
        
         vli::detail::single_coefficient_task empty_task;
         empty_task.step_count = 0;
         std::vector<vli::detail::single_coefficient_task> tasks_reordered(mul_block_size * max_iteration_count, empty_task);
    
         for(unsigned int batch_id = 0; batch_id < tasks.size() / 32; ++batch_id) {
                unsigned int warp_id = std::min_element(work_total_by_size.begin(), work_total_by_size.end()) - work_total_by_size.begin();
                std::copy(
                	tasks.begin() + (batch_id * 32),
                	tasks.begin() + ((batch_id + 1) * 32),
                	tasks_reordered.begin() + (workblock_count_by_warp_local[warp_id] * mul_block_size) + (warp_id * 32));
                unsigned int max_step_count = tasks[batch_id * 32].step_count;
        
                workblock_count_by_warp_local[warp_id]++;
                work_total_by_size[warp_id] += max_step_count;
         }
        
	 cudaMemcpyToSymbolAsync(workblock_count_by_warp, &(*workblock_count_by_warp_local.begin()), sizeof(unsigned int) * workblock_count_by_warp_local.size());
	 cudaMemcpyToSymbolAsync(execution_plan, &(*tasks_reordered.begin()), sizeof(vli::detail::single_coefficient_task) * tasks_reordered.size());
     }

     } // end namespace detail
 }//end namespace vli
