#include <boost/python.hpp>
#include "DmpAlgBgoDac.h"

BOOST_PYTHON_MODULE(libDmpAlgBgoDac){
  using namespace boost::python;

  class_<DmpAlgBgoDac,boost::noncopyable,bases<DmpVAlg> >("DmpAlgBgoDac",init<>());
}
