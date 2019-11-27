/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: classifyPoint_types.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 26-Nov-2019 23:34:49
 */

#ifndef CLASSIFYPOINT_TYPES_H
#define CLASSIFYPOINT_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_c_classreg_learning_coder_class
#define typedef_c_classreg_learning_coder_class

typedef struct {
  double CutVar[1003];
  double Children[2006];
  double ClassProb[29087];
  double CutPoint[1003];
  double ClassNames[29];
  int ClassNamesLength[29];
  double Prior[29];
  double Cost[841];
} c_classreg_learning_coder_class;

#endif                                 /*typedef_c_classreg_learning_coder_class*/
#endif

/*
 * File trailer for classifyPoint_types.h
 *
 * [EOF]
 */
