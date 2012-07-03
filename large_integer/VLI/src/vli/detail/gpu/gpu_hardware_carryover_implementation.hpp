
namespace vli {
    namespace detail {

    template <typename BaseInt, std::size_t Size, unsigned int Order>
    gpu_hardware_carryover_implementation<BaseInt, Size, Order>::gpu_hardware_carryover_implementation(){
        // As templated this array will be allocated a couple of time for every tupple of the cmake global size negligible  
        // only once due to singleton
        cudaMalloc((void**)&(this->execution_plan_), MulBlockSize<Order>::value*MaxIterationCount<Order>::value*sizeof(single_coefficient_task));
        cudaMalloc((void**)&(this->workblock_count_by_warp_), MulBlockSize<Order>::value/32*sizeof(unsigned int));
        element_count_prepared=0;
        plan();
    }

    template <typename BaseInt, std::size_t Size, unsigned int Order>
    gpu_hardware_carryover_implementation<BaseInt, Size, Order>::~gpu_hardware_carryover_implementation(){
        cudaFree(this->execution_plan_);
        cudaFree(this->workblock_count_by_warp_);
    } 

    template <typename BaseInt, std::size_t Size, unsigned int Order>
    void gpu_hardware_carryover_implementation<BaseInt, Size, Order>::plan(){
        std::vector<unsigned int> workblock_count_by_warp_local(MulBlockSize<Order>::value / 32,0);
        std::vector<unsigned int> work_total_by_size(MulBlockSize<Order>::value / 32,0);
        std::vector<vli::detail::single_coefficient_task> tasks((((Order*2+1)*(Order*2+1) + 32 - 1) / 32) * 32);
      
        for(unsigned int degree_y = 0; degree_y <(Order*2+1); ++degree_y) {
            for(unsigned int degree_x = 0; degree_x <(Order*2+1); ++degree_x) {
                vli::detail::single_coefficient_task& task = tasks[degree_y * (Order*2+1) + degree_x];
                task.output_degree_x = degree_x;
                task.output_degree_y = degree_y;
                task.step_count = (std::min<unsigned int>(((Order*2+1) - 1) - degree_x, degree_x) + 1) * (std::min<unsigned int>(((Order*2+1) - 1) - degree_y, degree_y) + 1);
            }
        }
        // Fill the task list up to the multiple of the warp size
        for(unsigned int i = (Order*2+1) * (Order*2+1); i < tasks.size(); ++i) {
               vli::detail::single_coefficient_task& task = tasks[i];
               task.output_degree_x = 0;
               task.output_degree_y = 0;
               task.step_count = 0;
        }
       // Sort the tasks in step_count descending order
         std::sort(tasks.begin(), tasks.end(), single_coefficient_task_sort);
        
         vli::detail::single_coefficient_task empty_task;
         empty_task.step_count = 0;
         std::vector<vli::detail::single_coefficient_task> tasks_reordered(MulBlockSize<Order>::value * MaxIterationCount<Order>::value, empty_task);
    
         for(unsigned int batch_id = 0; batch_id < tasks.size() / 32; ++batch_id) {
                unsigned int warp_id = std::min_element(work_total_by_size.begin(), work_total_by_size.end()) - work_total_by_size.begin();
                std::copy(
                	tasks.begin() + (batch_id * 32),
                	tasks.begin() + ((batch_id + 1) * 32),
                	tasks_reordered.begin() + (workblock_count_by_warp_local[warp_id] * MulBlockSize<Order>::value) + (warp_id * 32));
                unsigned int max_step_count = tasks[batch_id * 32].step_count;
        
                workblock_count_by_warp_local[warp_id]++;
                work_total_by_size[warp_id] += max_step_count;
         }
        
	 cudaMemcpyAsync(workblock_count_by_warp_, &(*workblock_count_by_warp_local.begin()), sizeof(unsigned int) * workblock_count_by_warp_local.size(), cudaMemcpyHostToDevice);
	 cudaMemcpyAsync(execution_plan_, &(*tasks_reordered.begin()), sizeof(single_coefficient_task) * tasks_reordered.size(),cudaMemcpyHostToDevice);
     }


     } // end namespace detail
 }//end namespace vli
