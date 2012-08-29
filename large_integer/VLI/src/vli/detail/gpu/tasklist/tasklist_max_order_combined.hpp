/*
 *Very Large Integer Library, License - Version 1.0 - May 3rd, 2012
 *
 *Timothee Ewart - University of Geneva,
 *Andreas Hehn - Swiss Federal Institute of technology Zurich.
 *Maxim Milakov – NVIDIA
 *
 *Permission is hereby granted, free of charge, to any person or organization
 *obtaining a copy of the software and accompanying documentation covered by
 *this license (the "Software") to use, reproduce, display, distribute,
 *execute, and transmit the Software, and to prepare derivative works of the
 *Software, and to permit third-parties to whom the Software is furnished to
 *do so, all subject to the following:
 *
 *The copyright notices in the Software and this entire statement, including
 *the above license grant, this restriction and the following disclaimer,
 *must be included in all copies of the Software, in whole or in part, and
 *all derivative works of the Software, unless such copies or derivative
 *works are solely in the form of machine-executable object code generated by
 *a source language processor.
 *
 *THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 *SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 *FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 *ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *DEALINGS IN THE SOFTWARE.
 */

#define VLI__ExtendStride  2*Order+1

namespace vli {
    namespace detail {
        // see doc for the origine of these equations
        // 2 variables            
        template< int Order>
        int CalculateStepCount(int output_degree_x,int output_degree_y){
             int sum(0);
            for(int i1 = std::max(0,(output_degree_x-Order)); i1 <= std::min(output_degree_x, Order); ++i1)
                sum += std::min(output_degree_y,(Order-i1)) - std::max(0,(output_degree_x+output_degree_y-Order-i1))+1;
            return sum;
            }        
       
        //3 variables
        template< int Order>
        int CalculateStepCount(int output_degree_x,int output_degree_y, int output_degree_z){
                 int sum(0);
                for(int i1 = std::max(0,(output_degree_x-Order)); i1 <= std::min(output_degree_x, Order); ++i1)
                    for(int i2 = std::max(0,(output_degree_x+output_degree_y-Order)-i1); i2 <= std::min(output_degree_y, Order-i1); ++i2)
                        sum += std::min(output_degree_z,(Order-i1-i2)) - std::max(0,(output_degree_x+output_degree_y+output_degree_z-Order-i1-i2))+1;
                return sum;
            }        

        //4 variables
        template< int Order>
        int CalculateStepCount(int output_degree_x,int output_degree_y, int output_degree_z, int output_degree_w){
                 int sum(0);
                for(int i1 = std::max(0,(output_degree_x-Order)); i1 <= std::min(output_degree_x, Order); ++i1)
                    for(int i2 = std::max(0,(output_degree_x+output_degree_y-Order)-i1); i2 <= std::min(output_degree_y, Order-i1); ++i2)
                        for(int i3 = std::max(0,(output_degree_x+output_degree_y+output_degree_z-Order)-i1-i2); i3 <= std::min(output_degree_z, Order-i1-i2); ++i3)
                            sum += std::min(output_degree_w,(Order-i1-i2-i3)) - std::max(0,(output_degree_x+output_degree_y+output_degree_z+output_degree_w-Order-i1-i2-i3))+1;
                return sum;
            }        
        
        template< int Order, int NumVars>
        struct BuildTaskList_helper;
        // one variable 'x', similar max_order_each one variable
        template< int Order>
        struct BuildTaskList_helper<Order,1>{
            static void BuildTaskList(std::vector<vli::detail::single_coefficient_task > & VecCoeff){
                for( int degree_x = 0; degree_x <VLI__ExtendStride; ++degree_x){
                    vli::detail::single_coefficient_task& task = VecCoeff[degree_x];
                    task.output_degree_x = degree_x;
                    task.output_degree_y = 0;
                    task.output_degree_z = 0;
                    task.output_degree_w = 0;
                    task.step_count =   (std::min< int>((VLI__ExtendStride - 1) - degree_x, degree_x) + 1);
                }
             }
        };
        // two variables 'x','y'
        template< int Order>
        struct BuildTaskList_helper<Order,2>{
            static void BuildTaskList(std::vector<vli::detail::single_coefficient_task > & VecCoeff){
                for( int degree_x = 0; degree_x <VLI__ExtendStride; ++degree_x)
                    for( int degree_y = 0; degree_y <VLI__ExtendStride - degree_x; ++degree_y){

                        vli::detail::single_coefficient_task& task = VecCoeff[VLI__ExtendStride*degree_x - (degree_x*degree_x-degree_x)/2 + degree_y];
                        task.output_degree_x = degree_x;
                        task.output_degree_y = degree_y;
                        task.output_degree_z = 0;
                        task.output_degree_w = 0;
                        task.step_count = CalculateStepCount<Order>(degree_x,degree_y);
                    }
            }
        };
        // tree variables 'x','y','z'
        template< int Order>
        struct BuildTaskList_helper<Order,3>{
            static void BuildTaskList(std::vector<vli::detail::single_coefficient_task > & VecCoeff){
                for( int degree_x = 0; degree_x <VLI__ExtendStride; ++degree_x)
                    for( int degree_y = 0; degree_y <VLI__ExtendStride - degree_x; ++degree_y)
                        for( int degree_z = 0; degree_z <VLI__ExtendStride - degree_x - degree_y; ++degree_z){

                            vli::detail::single_coefficient_task& task = VecCoeff[(degree_x*(degree_x*degree_x - 3*degree_x*(VLI__ExtendStride+1)
                                                                                  + 3*VLI__ExtendStride*(VLI__ExtendStride+2) +2))/6
                                                                                  + (VLI__ExtendStride - degree_x)*degree_y - (degree_y*degree_y-degree_y)/2 + degree_z];
                            task.output_degree_x = degree_x;
                            task.output_degree_y = degree_y;
                            task.output_degree_z = degree_z;
                            task.output_degree_w = 0;
                            task.step_count = CalculateStepCount<Order>(degree_x,degree_y,degree_z);
                    }
            }
        };
        // four variables 'x','y','z','w'
        template< int Order>
        struct BuildTaskList_helper<Order,4>{
            static void BuildTaskList(std::vector<vli::detail::single_coefficient_task > & VecCoeff){
                for( int degree_x = 0; degree_x <VLI__ExtendStride; ++degree_x)
                    for( int degree_y = 0; degree_y <VLI__ExtendStride - degree_x; ++degree_y)
                        for( int degree_z = 0; degree_z <VLI__ExtendStride - degree_x - degree_y; ++degree_z)
                            for( int degree_w = 0; degree_w <VLI__ExtendStride - degree_x - degree_y - degree_z; ++degree_w){

                                vli::detail::single_coefficient_task& task = VecCoeff[(degree_x*(2*VLI__ExtendStride+3-degree_x)*(2*VLI__ExtendStride*VLI__ExtendStride+6*VLI__ExtendStride+2 +degree_x*degree_x -2*VLI__ExtendStride*degree_x - 3*degree_x))/24
                                                                                     + (degree_y*(degree_y*degree_y - 3*degree_y*(VLI__ExtendStride+1-degree_x) + 3*(VLI__ExtendStride-degree_x)*(VLI__ExtendStride+2-degree_x)+2))/6
                                                                                     +(VLI__ExtendStride-degree_x-degree_y)*degree_z - (degree_z*degree_z-degree_z)/2 + degree_w];
                                task.output_degree_x = degree_x;
                                task.output_degree_y = degree_y;
                                task.output_degree_z = degree_z;
                                task.output_degree_w = degree_w;
                                task.step_count = CalculateStepCount<Order>(degree_x,degree_y,degree_z,degree_w);

                            }
            }
        };

    template <std::size_t NumBits, int Order, int NumVars>
    tasklist_keep_order<NumBits, max_order_combined<Order>, NumVars>::tasklist_keep_order(){
        // As templated this array will be allocated a couple of time for every tupple of the cmake global size negligible  
        // only once due to singleton
        gpu::cu_check_error(cudaMalloc((void**)&(this->execution_plan_), mul_block_size<max_order_combined<Order>, NumVars>::value*MaxIterationCount<max_order_combined<Order>, NumVars>::value*sizeof(single_coefficient_task)),__LINE__);
        gpu::cu_check_error(cudaMalloc((void**)&(this->workblock_count_by_warp_), mul_block_size<max_order_combined<Order>, NumVars>::value/32*sizeof( int)),__LINE__);
        element_count_prepared=0;
        plan();
    }

    template <std::size_t NumBits, int Order, int NumVars>
    tasklist_keep_order<NumBits, max_order_combined<Order>, NumVars>::~tasklist_keep_order(){
        gpu::cu_check_error(cudaFree(this->execution_plan_),__LINE__);
        gpu::cu_check_error(cudaFree(this->workblock_count_by_warp_),__LINE__);
    } 

    template <std::size_t NumBits, int Order, int NumVars>
    void tasklist_keep_order<NumBits, max_order_combined<Order>, NumVars>::plan(){
        std::vector<int> workblock_count_by_warp_local(mul_block_size<max_order_combined<Order>, NumVars>::value / 32U,0);
        std::vector<int> work_total_by_size(mul_block_size<max_order_combined<Order>, NumVars>::value / 32U,0);
        //TO CHECK 2*Order or 2*order+1
        std::vector<vli::detail::single_coefficient_task > tasks(((vli::detail::max_order_combined_helpers::size<NumVars+1, 2*Order>::value
                                                                   + 32U - 1) / 32U) * 32U);
        
        BuildTaskList_helper<Order,NumVars>::BuildTaskList(tasks);
        
        // Fill the task list up to the multiple of the warp size
        for(unsigned  int i = vli::detail::max_order_combined_helpers::size<NumVars+1, 2*Order>::value; i < tasks.size(); ++i) {
               vli::detail::single_coefficient_task& task = tasks[i];
               task.output_degree_x = 0;
               task.output_degree_y = 0;
               task.output_degree_z = 0;
               task.output_degree_w = 0;
               task.step_count = 0;
        }
       // Sort the tasks in step_count descending order
         std::sort(tasks.begin(), tasks.end(), vli::detail::single_coefficient_task_sort);
         std::vector<vli::detail::single_coefficient_task > tasks_reordered(mul_block_size<max_order_combined<Order>, NumVars>::value
                                                                            * MaxIterationCount<max_order_combined<Order>, NumVars>::value);
         // this thing should be generic ... yes it is ! 
         for( unsigned int batch_id = 0; batch_id < tasks.size() / 32; ++batch_id) {
                 //TO DO : std::distance more safe !!!!!
                 int warp_id = std::min_element(work_total_by_size.begin(), work_total_by_size.end()) - work_total_by_size.begin(); // - to get the position
                std::copy(
                	tasks.begin() + (batch_id * 32),
                	tasks.begin() + ((batch_id + 1) * 32),
                	tasks_reordered.begin() + (workblock_count_by_warp_local[warp_id] * mul_block_size<max_order_combined<Order>, NumVars>::value) + (warp_id * 32));
                 int max_step_count = tasks[batch_id * 32].step_count;
        
                workblock_count_by_warp_local[warp_id]++;
                work_total_by_size[warp_id] += max_step_count;
         }
	 gpu::cu_check_error(cudaMemcpyAsync(workblock_count_by_warp_, &(*workblock_count_by_warp_local.begin()), sizeof( int) * workblock_count_by_warp_local.size(), cudaMemcpyHostToDevice),__LINE__);
         gpu::cu_check_error(cudaMemcpyAsync(execution_plan_, &(*tasks_reordered.begin()), sizeof(single_coefficient_task) * tasks_reordered.size(),cudaMemcpyHostToDevice),__LINE__);
    }

    } // end namespace detail
 }//end namespace vli
#undef VLI__ExtendStride // delete the macro
