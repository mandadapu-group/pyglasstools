#ifndef __LOGFILE_H__
#define __LOGFILE_H__

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>

#include <pybind11/pybind11.h>


/* Class: BaseLogFile
 * Provides the base class for read and write of logfiles
 */
class PYBIND11_EXPORT BaseLogFile
{
    public:
        BaseLogFile()
        {
        };
        
        //! Destructor
        virtual ~BaseLogFile()
        {
        };

        virtual void write(std::string input)
        {
        };
};



/* Class: LogFile
 * The class for read and write of logfiles
 * when MPI is not used
 */
class PYBIND11_EXPORT LogFile : public BaseLogFile
{
    private:
        std::fstream m_file;
    public:
        LogFile(std::string filename)
        {
            m_file.open(filename, std::ios::out | std::ios::app);
        };
        
        //! Destructor
        virtual ~LogFile()
        {
            m_file.close();
        };

        void write(std::string input)
        {
            m_file << input;
        };
};

//! Exports Communicator to python
void export_BaseLogFile(pybind11::module& m)
{
    pybind11::class_< BaseLogFile, std::shared_ptr<BaseLogFile> >(m,"BaseLogFile")
    .def("write", &BaseLogFile::write)
    ;
};

//! Exports Communicator to python
void export_LogFile(pybind11::module& m)
{
    pybind11::class_< LogFile, BaseLogFile, std::shared_ptr<LogFile> >(m,"LogFile")
    .def( pybind11::init< std::string >() )
    .def("write", &LogFile::write)
    ;
};

#endif
