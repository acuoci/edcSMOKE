/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: classifyPoint.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 26-Nov-2019 23:34:49
 */

/* Include Files */
#include <string.h>
#include "rt_nonfinite.h"
#include "classifyPoint.h"
#include "CompactClassificationTree.h"
#include "CompactClassificationModel.h"
#include "CompactTree.h"

/* Function Definitions */

/*
 * Arguments    : const double X[6]
 * Return Type  : double
 */
double classifyPoint(const double X[6])
{
  static c_classreg_learning_coder_class CompactMdl;
  static const double dv0[29] = { 0.14913, 0.01561, 0.03713, 0.0006, 9.5E-5,
    0.042275, 0.006845, 0.012685, 0.012025, 6.0E-5, 0.00793, 0.00045, 0.000755,
    0.00615, 0.049595, 0.000835, 0.005815, 6.5E-5, 0.015915, 0.048955, 0.012865,
    0.001875, 0.00546, 0.004665, 0.15609, 1.0E-5, 0.002765, 0.03943, 0.36392 };

  c_CompactClassificationModel_Co(&CompactMdl);
  CompactTree_CompactTree(&CompactMdl);
  memcpy(&CompactMdl.Prior[0], &dv0[0], 29U * sizeof(double));
  c_CompactClassificationModel_se(&CompactMdl);
  return c_CompactClassificationTree_pre(CompactMdl.CutVar, CompactMdl.Children,
    CompactMdl.ClassProb, CompactMdl.CutPoint, CompactMdl.ClassNames,
    CompactMdl.Cost, X);
}

/*
 * File trailer for classifyPoint.c
 *
 * [EOF]
 */
