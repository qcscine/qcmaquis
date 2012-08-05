#define VLI__ExtendStride extend_stride<Var0, Order>::value // 2*order+1

namespace vli {
    namespace detail {
        
        
        template< unsigned int Order, class Var0, class Var1, class Var2, class Var3>
        struct BuildTaskList_helper;
        // one variable 'x', similar max_order_each one variable
        template< unsigned int Order, class Var0>
        struct BuildTaskList_helper<Order, Var0, vli::no_variable, vli::no_variable, vli::no_variable>{
            static void BuildTaskList(std::vector<vli::detail::single_coefficient_task > & VecCoeff){
                for(unsigned int degree_x = 0; degree_x <VLI__ExtendStride; ++degree_x){
                    vli::detail::single_coefficient_task& task = VecCoeff[degree_x];
                    task.output_degree_x = degree_x;
                    task.output_degree_y = 0;
                    task.output_degree_z = 0;
                    task.output_degree_w = 0;
                    task.step_count =   (std::min<unsigned int>((VLI__ExtendStride - 1) - degree_x, degree_x) + 1);
                }
             }
        };
        // two variables 'x','y'
        template< unsigned int Order, class Var0, class Var1>
        struct BuildTaskList_helper<Order, Var0, Var1, vli::no_variable, vli::no_variable>{
            static void BuildTaskList(std::vector<vli::detail::single_coefficient_task > & VecCoeff){
                for(unsigned int degree_x = 0; degree_x <VLI__ExtendStride; ++degree_x)
                    for(unsigned int degree_y = 0; degree_y <VLI__ExtendStride - degree_x; ++degree_y){
                        //TO CHECK 2*order, come from Andreas work
                        vli::detail::single_coefficient_task& task = VecCoeff[VLI__ExtendStride*degree_x - (degree_x*degree_x-degree_x)/2 + degree_y];
                        task.output_degree_x = degree_x;
                        task.output_degree_y = degree_y;
                        task.output_degree_z = 0;
                        task.output_degree_w = 0;
                        //TO CHECK is it thrue ?
                        task.step_count =   (std::min<unsigned int>((VLI__ExtendStride - 1) - degree_x, degree_x) + 1) // C EST TOTALEMENT FAUX !@!!!!!!!!!asdjhvbaksdhgqvksdvbqoy, WRONG TO CHANGE !!!
                                          * (std::min<unsigned int>((VLI__ExtendStride - 1) - degree_y, degree_y) + 1) 
                    }
            }
        };
        // tree variables 'x','y','z'
        template< unsigned int Order, class Var0, class Var1, class Var2>
        struct BuildTaskList_helper<Order, Var0, Var1, Var2, vli::no_variable>{
            static void BuildTaskList(std::vector<vli::detail::single_coefficient_task > & VecCoeff){
                for(unsigned int degree_x = 0; degree_x <VLI__ExtendStride; ++degree_x)
                    for(unsigned int degree_y = 0; degree_y <VLI__ExtendStride - degree_x; ++degree_y)
                        for(unsigned int degree_z = 0; degree_z <VLI__ExtendStride - degree_x - degree_y; ++degree_z){
                            //TO CHECK 2*order, come from Andreas work
                            vli::detail::single_coefficient_task& task = VecCoeff[(degree_x*(degree_x*degree_x - 3*degree_x*(VLI__ExtendStride+1)
                                                                                  + 3*VLI__ExtendStride*(VLI__ExtendStride+2) +2))/6
                                                                                  + (VLI__ExtendStride - degree_x)*degree_y - (degree_y*degree_y-degree_y)/2 + degree_z];
                            task.output_degree_x = degree_x;
                            task.output_degree_y = degree_y;
                            task.output_degree_z = degree_z;
                            task.output_degree_w = 0;
                            //TO CHECK is it thrue ?
                            task.step_count =   (std::min<unsigned int>((VLI__ExtendStride - 1) - degree_x, degree_x) + 1)
                                              * (std::min<unsigned int>((VLI__ExtendStride - 1) - degree_y, degree_y) + 1) 
                                              * (std::min<unsigned int>((VLI__ExtendStride - 1) - degree_z, degree_z) + 1) ;
                    }
            }
        };
        // four variables 'x','y','z','w'
        template< unsigned int Order, class Var0, class Var1, class Var2, class Var3>
        struct BuildTaskList_helper{
            static void BuildTaskList(std::vector<vli::detail::single_coefficient_task > & VecCoeff){
                for(unsigned int degree_x = 0; degree_x <VLI__ExtendStride; ++degree_x)
                    for(unsigned int degree_y = 0; degree_y <VLI__ExtendStride - degree_x; ++degree_y)
                        for(unsigned int degree_z = 0; degree_z <VLI__ExtendStride - degree_x - degree_y; ++degree_z)
                            for(unsigned int degree_w = 0; degree_w <VLI__ExtendStride - degree_x - degree_y - degree_y; ++degree_w){
                                //TO CHECK 2*order, come from Andreas work
                                vli::detail::single_coefficient_task& task = VecCoeff[(degree_x*(2*VLI__ExtendStride+3-degree_x)*(2*VLI__ExtendStride*VLI__ExtendStride+6*VLI__ExtendStride+2 +degree_x*degree_x -2*VLI__ExtendStride*degree_x - 3*degree_x))/24
                                                                                      + (degree_y*(degree_y*degree_y - 3*degree_y*(VLI__ExtendStride+1-degree_x) + 3*(VLI__ExtendStride-degree_x)*(VLI__ExtendStride+2-degree_x)+2))/6
                                                                                      +(VLI__ExtendStride-degree_x-degree_y)*degree_z - (degree_z*degree_z-degree_z)/2 + degree_w];
                                                                                      
                                task.output_degree_x = degree_x;
                                task.output_degree_y = degree_y;
                                task.output_degree_z = degree_z;
                                task.output_degree_w = degree_w;
                                //TO CHECK is it thrue ?
                                task.step_count =   (std::min<unsigned int>((extend_stride<Var0, Order>::value - 1) - degree_x, degree_x) + 1)
                                                  * (std::min<unsigned int>((extend_stride<Var1, Order>::value - 1) - degree_y, degree_y) + 1) 
                                                  * (std::min<unsigned int>((extend_stride<Var2, Order>::value - 1) - degree_z, degree_z) + 1)
                                                  * (std::min<unsigned int>((extend_stride<Var3, Order>::value - 1) - degree_w, degree_w) + 1) ;
                            }
            }
        };

    template <std::size_t Size, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    tasklist_keep_order<Size, max_order_combined<Order>, Var0, Var1, Var2, Var3>::tasklist_keep_order(){
        // As templated this array will be allocated a couple of time for every tupple of the cmake global size negligible  
        // only once due to singleton
        cudaMalloc((void**)&(this->execution_plan_), MulBlockSize<max_order_combined<Order>, Var0, Var1, Var2, Var3>::value*MaxIterationCount<max_order_combined<Order>, Var0, Var1, Var2, Var3>::value*sizeof(single_coefficient_task));
        cudaMalloc((void**)&(this->workblock_count_by_warp_), MulBlockSize<max_order_combined<Order>, Var0, Var1, Var2, Var3>::value/32*sizeof(unsigned int));
        element_count_prepared=0;
        plan();
    }

    template <std::size_t Size, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    tasklist_keep_order<Size, max_order_combined<Order>, Var0, Var1, Var2, Var3>::~tasklist_keep_order(){
        cudaFree(this->execution_plan_);
        cudaFree(this->workblock_count_by_warp_);
    } 

    template <std::size_t Size, unsigned int Order, class Var0, class Var1, class Var2, class Var3>
    void tasklist_keep_order<Size, max_order_combined<Order>, Var0, Var1, Var2, Var3>::plan(){
        std::vector<unsigned int> workblock_count_by_warp_local(MulBlockSize<max_order_combined<Order>, Var0, Var1, Var2, Var3>::value / 32U,0);
        std::vector<unsigned int> work_total_by_size(MulBlockSize<max_order_combined<Order>, Var0, Var1, Var2, Var3>::value / 32U,0);
        //TO CHECK 2*Order or 2*order+1
        std::vector<vli::detail::single_coefficient_task > tasks(((vli::detail::max_order_combined_helpers::size<vli::detail::num_of_variables_helper<Var0,Var1,Var2,Var3 >::value+1, 2*Order>::value
                                                                   + 32U - 1) / 32U) * 32U);
        
        BuildTaskList_helper<Order,Var0,Var1,Var2,Var3>::BuildTaskList(tasks);
        
        // Fill the task list up to the multiple of the warp size
        for(unsigned int i = vli::detail::max_order_combined_helpers::size<vli::detail::num_of_variables_helper<Var0,Var1,Var2,Var3 >::value+1, 2*Order>::value; i < tasks.size(); ++i) {
               vli::detail::single_coefficient_task& task = tasks[i];
               task.output_degree_x = 0;
               task.output_degree_y = 0;
               task.output_degree_z = 0;
               task.output_degree_w = 0;
               task.step_count = 0;
        }
       // Sort the tasks in step_count descending order
         std::sort(tasks.begin(), tasks.end(), vli::detail::single_coefficient_task_sort);
         std::vector<vli::detail::single_coefficient_task > tasks_reordered(MulBlockSize<max_order_combined<Order>, Var0, Var1, Var2, Var3>::value
                                                                            * MaxIterationCount<max_order_combined<Order>, Var0, Var1, Var2, Var3>::value);
         // this thing should be generic ... yes it is ! 
         for(unsigned int batch_id = 0; batch_id < tasks.size() / 32; ++batch_id) {
                unsigned int warp_id = std::min_element(work_total_by_size.begin(), work_total_by_size.end()) - work_total_by_size.begin(); // - to get the position
                std::copy(
                	tasks.begin() + (batch_id * 32),
                	tasks.begin() + ((batch_id + 1) * 32),
                	tasks_reordered.begin() + (workblock_count_by_warp_local[warp_id] * MulBlockSize<max_order_combined<Order>, Var0, Var1, Var2, Var3>::value) + (warp_id * 32));
                unsigned int max_step_count = tasks[batch_id * 32].step_count;
        
                workblock_count_by_warp_local[warp_id]++;
                work_total_by_size[warp_id] += max_step_count;
         }
       
	 cudaMemcpyAsync(workblock_count_by_warp_, &(*workblock_count_by_warp_local.begin()), sizeof(unsigned int) * workblock_count_by_warp_local.size(), cudaMemcpyHostToDevice);
         cudaMemcpyAsync(execution_plan_, &(*tasks_reordered.begin()), sizeof(single_coefficient_task) * tasks_reordered.size(),cudaMemcpyHostToDevice);
    }

    } // end namespace detail
 }//end namespace vli
#undef VLI__ExtendStride // delete the macro
