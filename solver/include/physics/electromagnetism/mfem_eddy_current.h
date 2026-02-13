#pragma once

#include "config.h"

#ifdef LOAD_MFEM
  #include "mfem.hpp"
#endif

class MFEM_Eddy_Current
{

private:
  const char * mesh_file_;
  int order_;
  bool pa_;
  const char * device_config_;
  bool visualization_;

public:

  MFEM_Eddy_Current(const char * mesh_file, int order=1, bool pa=false, const char * device_config="cpu", bool visualization=1);
};