/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: CompactClassificationModel.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 26-Nov-2019 23:34:49
 */

/* Include Files */
#include <string.h>
#include "rt_nonfinite.h"
#include "classifyPoint.h"
#include "CompactClassificationModel.h"

/* Function Definitions */

/*
 * Arguments    : c_classreg_learning_coder_class *obj
 * Return Type  : void
 */
void c_CompactClassificationModel_Co(c_classreg_learning_coder_class *obj)
{
  int i;
  for (i = 0; i < 29; i++) {
    obj->ClassNames[i] = 1.0 + (double)i;
    obj->ClassNamesLength[i] = 1;
  }
}

/*
 * Arguments    : c_classreg_learning_coder_class *obj
 * Return Type  : void
 */
void c_CompactClassificationModel_se(c_classreg_learning_coder_class *obj)
{
  signed char I[841];
  int k;
  memset(&I[0], 0, 841U * sizeof(signed char));
  for (k = 0; k < 29; k++) {
    I[k + 29 * k] = 1;
  }

  for (k = 0; k < 841; k++) {
    obj->Cost[k] = 1.0 - (double)I[k];
  }
}

/*
 * File trailer for CompactClassificationModel.c
 *
 * [EOF]
 */
