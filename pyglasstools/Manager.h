#ifndef __MANAGER_H__
#define __MANAGER_H__

#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <mpi.h>

#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/embed.h>
#include <float.h>
namespace py = pybind11;

namespace detail
{
    template <typename T>
    void argparser(const char command[], T &val, std::string errormsg, std::map< std::string, std::string > argv);

    template <>
    void argparser<double>(const char command[], double &val, std::string errormsg, std::map< std::string, std::string > argv)
    {
        try
        {
            if (!argv[command].empty())
            {
                val = std::stod(argv[std::string(command)],nullptr);
            }
        }
        catch(...)
        {
            throw std::domain_error(errormsg);
        }
    };

    template <>
    void argparser<int>(const char command[], int &val, std::string errormsg, std::map< std::string, std::string > argv)
    {
        try
        {
            if (!argv[command].empty())
            {
                val = std::stoi(argv[command],nullptr);
            }
        }
        catch(...)
        {
            throw std::domain_error(errormsg);
        }
    };
    
    template <>
    void argparser<unsigned int>(const char command[], unsigned int &val, std::string errormsg, std::map< std::string, std::string > argv)
    {
        try
        {
            if (!argv[command].empty())
            {
                val = (unsigned int)std::stoi(argv[command],nullptr);
            }
        }
        catch(...)
        {
            throw std::domain_error(errormsg);
        }
    };

    template <>
    void argparser<std::string>(const char command[], std::string &val, std::string errormsg, std::map< std::string, std::string > argv)
    {

        try
        {
            if (!argv[command].empty())
            {
                val = argv[command];
            }
        }
        catch(...)
        {
            throw std::domain_error(errormsg);
        }
    };

    void get_argv(std::map< std::string, std::string >& argv_list, std::string& cmd_line_options)
    {
        auto argv = py::module::import("sys").attr("argv");
        for( auto test = argv.begin(); test != argv.end(); ++test)
        {
            std::string s1 = std::string(py::str(*test));
            //specific to petsc/slepc related commands
            if (s1.at(0) == '-')
            {
                std::string val = std::string(py::str(*std::next(test,1)));
                //if (val.at(0) != '-')
                //{
                argv_list.insert(std::pair<std::string, std::string >(s1, val));
                cmd_line_options += s1+" "+val+" ";
                //}
            }
        }
    }

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
};


class PYBIND11_EXPORT Manager
{
    protected:
        std::map<std::string, std::string> argv;
    public:
        //Basic MPI information
        std::string cmd_line_options;  
        int proc_rank, nprocs;
        
        Manager()  
        {
            m_notice_prefix  = "notice";
            m_nullstream = std::shared_ptr< detail::nullstream >(new detail::nullstream());
            m_notice_stream = &std::cout;
            MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

            //Default value of notice_level is 1
            m_notice_level = 1;        
            detail::get_argv(argv,cmd_line_options);//py::module::import("sys").attr("argv");
            detail::argparser<unsigned int>("-notice_level",m_notice_level, "[ERROR] Invalid value for notice_level", argv);
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
        std::ostream& notice(unsigned int level, int target_rank = 0)
        {
            assert(m_notice_stream);
            if (level <= m_notice_level && proc_rank == target_rank)
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
        void py_notice(unsigned int level, std::string message)
        {
            std::stringstream notice_stream;//assert(m_notice_stream);
            if (level <= m_notice_level && proc_rank == 0 && m_notice_prefix != std::string("") && level > 1)
            {
                notice_stream << m_notice_prefix << "(" << level << "): " << message;
                py::print(notice_stream.str());
            }
        }
        std::ostream& widenotice(unsigned int level)
        {
            assert(m_notice_stream);
            if (level <= m_notice_level)// && proc_rank == 0)
                {
                if (m_notice_prefix != std::string("") && level > 1)
                    *m_notice_stream << m_notice_prefix << "(" << level << ") for process (" << proc_rank << "): ";
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
        std::shared_ptr< detail::nullstream >  m_nullstream;   //!< null stream
        unsigned int m_notice_level;
};

void export_Manager(py::module& m)
{
    py::class_<Manager, std::shared_ptr<Manager> >(m,"Manager")
    .def(py::init<>())
    .def("notice",&Manager::py_notice)
    ;
};
#endif
