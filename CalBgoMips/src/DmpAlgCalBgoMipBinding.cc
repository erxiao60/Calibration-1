/*
 *  $Id: DmpAlgCalBgoMipBinding.cc, 2014-09-03 11:35:31 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 03/09/2014
*/

#include <boost/python.hpp>
#include "DmpAlgBgoMips.h"

BOOST_PYTHON_MODULE(libDmpBgoMips){
  using namespace boost::python;

  class_<DmpAlgBgoMips,boost::noncopyable,bases<DmpVAlg> >("DmpAlgBgoMips",init<>());
}

