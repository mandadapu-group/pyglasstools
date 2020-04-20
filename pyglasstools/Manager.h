#ifndef __MANAGER_H__
#define __MANAGER_H__

#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <mpi.h>

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/embed.h>
#include <float.h>
namespace py = pybind11;

//! A null stream that doesn't write anything sent to it
/*! From: http://bytes.com/topic/c/answers/127843-null-output-stream#post444998
*/
struct PYBIND11_EXPORT nullstream: std::ostream
    {
    //! Construct a null stream
    nullstream(): std::ios(0), std::ostream(0) {}
    };

template <typename T>
std::string to_string_sci(T value) 
{
      std::stringstream out;
      //out << std::scientific;
      out << std::setprecision(std::numeric_limits<double>::max_digits10-1);
      out << value;
    return out.str();
}

class PYBIND11_EXPORT Manager
{
    public:
        //Basic MPI information
        int proc_rank, nprocs;
        std::string cmd_line_options; 
        
        Manager()  
        {
            m_notice_prefix  = "notice";
            m_nullstream = std::shared_ptr<nullstream>(new nullstream());
            m_notice_stream = &std::cout;
            MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

            //Default value of notice_level is 1
            m_notice_level = 1;        
            auto argv = py::module::import("sys").attr("argv");
            int argc = 0;
            for( auto test = argv.begin(); test != argv.end(); ++test)
            {
                if (std::string(py::str(*test)) == "-notice_level")
                {
                    //Increment manually
                    ++test;
                    try
                    {
                        m_notice_level = std::stoi(py::str(*test),nullptr);
                    }
                    catch(...)
                    {
                        throw std::domain_error(std::string("Process [")+std::to_string(proc_rank)+std::string("]: ") +
                                                std::string("Did you correctly input your notice level? The value you put is: ")+
                                                std::string(py::str(*test)));
                    }
                }
                else
                {
                    if (argc == 0)
                        ++argc;
                    else
                    {
                        cmd_line_options += py::str(*test);
                        cmd_line_options += " ";
                        argc += 1;
                    }
                }
            }
        };
        virtual ~Manager()
        {
            m_notice_stream = NULL;
        };  
        void setNoticeLevel(int notice_level)
        {
            m_notice_level = notice_level;
        } 
        unsigned int getNoticeLevel()
        {
            return m_notice_level;
        } 
        std::ostream& getStream()
        {
            return *m_notice_stream;
        } 
        std::ostream& notice(unsigned int level)
        {
            assert(m_notice_stream);
            if (level <= m_notice_level && proc_rank == 0)
                {
                if (m_notice_prefix != std::string("") && level > 1)
                    *m_notice_stream << m_notice_prefix << "(" << level << "): ";
                return *m_notice_stream;
                }
            else
                {
                return *m_nullstream;
                }
        }
    protected:
        std::string m_notice_prefix;
        std::ostream *m_notice_stream;  //!< notice stream
        std::shared_ptr<nullstream>  m_nullstream;   //!< null stream
        unsigned int m_notice_level;
};

#endif
