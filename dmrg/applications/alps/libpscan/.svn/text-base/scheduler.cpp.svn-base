/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "libpscan/scheduler.hpp"
#include "libpscan/run_sim.hpp"

#include "dmrg/utils/time_limit_exception.h"

#include <alps/parser/parser.h>
#include <alps/parser/xslt_path.h>
#include <alps/parser/xmlstream.h>

#include <boost/filesystem/fstream.hpp>
#include <iostream>

Scheduler::Scheduler(const Options& opt)
: stop_callback(opt.time_limit)
{
    outfilepath = opt.jobfilename;
    infilepath = opt.jobfilename;
    
    infilepath = boost::filesystem::absolute(infilepath);
    outfilepath = boost::filesystem::absolute(outfilepath);
    
    if (!opt.jobfilename.empty())
        parse_job_file(infilepath);
}

void Scheduler::parse_job_file(const boost::filesystem::path& filename)
{
    boost::filesystem::ifstream infile(filename);
    alps::XMLTag tag = alps::parse_tag(infile, true);
    if (tag.name != "JOB")
        boost::throw_exception(std::runtime_error("missing <JOB> element in jobfile"));
    tag = alps::parse_tag(infile);
    if (tag.name == "OUTPUT") {
        if(tag.attributes["file"] != "")
            outfilepath = boost::filesystem::absolute(boost::filesystem::path(tag.attributes["file"]), filename.parent_path());
        else
            boost::throw_exception(std::runtime_error("missing 'file' attribute in <OUTPUT> element in jobfile"));
        tag = alps::parse_tag(infile);
        if (tag.name == "/OUTPUT")
            tag = alps::parse_tag(infile);
    }
    
    while (tag.name == "TASK") {
        TaskDescriptor task;
        if (tag.attributes["status"] == "" || tag.attributes["status"] == "new")
            task.status = TaskNotStarted;
        else if (tag.attributes["status"] == "running")
            task.status = TaskHalted;
        else if (tag.attributes["status"] == "finished")
            task.status = TaskFinished;
        else
            boost::throw_exception(std::runtime_error("illegal status attribute in <TASK> element in jobfile"));
        
        tag = alps::parse_tag(infile);
        if (tag.name == "INPUT") {
            task.in = boost::filesystem::path(tag.attributes["file"]);
            if (task.in.empty())
                boost::throw_exception(std::runtime_error("missing 'file' attribute in <INPUT> element in jobfile"));
            tag = alps::parse_tag(infile);
            if (tag.name == "/INPUT")
                tag = alps::parse_tag(infile);
        }
        else
            boost::throw_exception(std::runtime_error("missing <INPUT> element in jobfile"));
        
        if (tag.name == "OUTPUT") {
            task.out = boost::filesystem::path(tag.attributes["file"]);
            if (task.out.empty())
                boost::throw_exception(std::runtime_error("missing 'file' attribute in <OUTPUT> element in jobfile"));
            tag = alps::parse_tag(infile);
            if (tag.name == "/OUTPUT")
                tag = alps::parse_tag(infile);
        }
        if (task.out.empty())
            task.out = task.in;
        task.in  = boost::filesystem::absolute(task.in, filename.parent_path());
        task.out = boost::filesystem::absolute(task.out, filename.parent_path());
        if (tag.name != "/TASK")
            boost::throw_exception(std::runtime_error("missing </TASK> tag in jobfile"));
        tag = alps::parse_tag(infile);
        tasks.push_back(task);
    }
    if (tag.name != "/JOB")
        boost::throw_exception(std::runtime_error("missing </JOB> tag in jobfile"));
}


void Scheduler::checkpoint_status() const
{
    bool make_backup = boost::filesystem::exists(outfilepath);
    boost::filesystem::path filename = outfilepath;
    boost::filesystem::path dir = outfilepath.parent_path();
    if (make_backup)
        filename = dir / (filename.filename().string() + ".bak");
    { // scope for out
        alps::oxstream out(filename);
        
        out << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("ALPS.xsl"));
        out << alps::start_tag("JOB")
            << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
            << alps::attribute("xsi:noNamespaceSchemaLocation",
                               "http://xml.comp-phys.org/2003/8/job.xsd");
        
        for (unsigned int i=0; i<tasks.size(); i++) {
            std::string status;
            if (tasks[i].status == TaskFinished)
                status = "finished";
            else if (tasks[i].status == TaskNotStarted)
                status = "new";
            else // TaskRunning or TaskHalted
                status = "running";
            
            if (tasks[i].status == TaskNotStarted && !boost::filesystem::exists(tasks[i].out))
                copy(tasks[i].in, tasks[i].out);
            
            out << alps::start_tag("TASK") << alps::attribute("status", status)
                << alps::start_tag("INPUT")
                << alps::attribute("file", tasks[i].out.string())
                << alps::end_tag() << alps::end_tag();
        }
        out << alps::end_tag("JOB");
    }
    if (make_backup) {
        boost::filesystem::remove(outfilepath);
        boost::filesystem::rename(filename, outfilepath);
    }
}


void Scheduler::run()
{
    // do all Tasks
    try {
        for(unsigned int i=0; i<tasks.size(); i++) {
            boost::chrono::high_resolution_clock::time_point t0 = boost::chrono::high_resolution_clock::now();
            double time_left = stop_callback.time_left().count();
            if (stop_callback.valid() && time_left < 0)
                throw dmrg::time_limit();
                
            if (tasks[i].status == TaskFinished)
                std::cout << "Task " << i+1 << " finished." << std::endl;
            else if (tasks[i].status == TaskNotStarted || tasks[i].status == TaskRunning ||
                     tasks[i].status == TaskHalted) {
                
                tasks[i].status = TaskRunning;
                std::cout  << "Running task " << i+1 << "." << std::endl;
                if (!boost::filesystem::exists(tasks[i].out))
                    copy(tasks[i].in, tasks[i].out);

                /// Start task
                run_sim(tasks[i].in, tasks[i].out, time_left);
                tasks[i].status = TaskFinished;
            }
            else
                boost::throw_exception( std::logic_error("illegal Task status"));
            
            boost::chrono::high_resolution_clock::time_point t1 = boost::chrono::high_resolution_clock::now();
            std::cout  << "This task took " << boost::chrono::duration<double>(t1-t0) << "." << std::endl;
            checkpoint_status();
        }
        std::cout << "Finished with everything." << std::endl;
    } catch (dmrg::time_limit const& e) {
        std::cout << "Time limit exceeded." << std::endl;
        checkpoint_status();
    }
}
