/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: CompactClassificationTree.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 26-Nov-2019 23:34:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "classifyPoint.h"
#include "CompactClassificationTree.h"

/* Function Definitions */

/*
 * Arguments    : const double obj_CutVar[1003]
 *                const double obj_Children[2006]
 *                const double obj_ClassProb[29087]
 *                const double obj_CutPoint[1003]
 *                const double obj_ClassNames[29]
 *                const double obj_Cost[841]
 *                const double X[6]
 * Return Type  : double
 */
double c_CompactClassificationTree_pre(const double obj_CutVar[1003], const
  double obj_Children[2006], const double obj_ClassProb[29087], const double
  obj_CutPoint[1003], const double obj_ClassNames[29], const double obj_Cost[841],
  const double X[6])
{
  double m;
  static const unsigned char pruneList[1003] = { 139U, 138U, 139U, 136U, 121U,
    133U, 137U, 132U, 113U, 121U, 42U, 120U, 62U, 116U, 135U, 127U, 130U, 113U,
    113U, 112U, 121U, 7U, 42U, 0U, 87U, 26U, 31U, 73U, 116U, 118U, 135U, 115U,
    99U, 130U, 102U, 101U, 108U, 75U, 113U, 94U, 71U, 104U, 121U, 7U, 7U, 0U, 0U,
    18U, 60U, 0U, 0U, 28U, 0U, 73U, 73U, 95U, 85U, 76U, 118U, 131U, 134U, 50U,
    115U, 38U, 93U, 126U, 124U, 0U, 98U, 65U, 56U, 39U, 73U, 39U, 75U, 108U, 10U,
    36U, 39U, 40U, 71U, 104U, 51U, 23U, 105U, 0U, 7U, 0U, 0U, 0U, 0U, 53U, 60U,
    28U, 13U, 0U, 32U, 73U, 0U, 82U, 21U, 0U, 85U, 35U, 76U, 61U, 53U, 128U,
    129U, 122U, 117U, 47U, 0U, 83U, 96U, 38U, 14U, 0U, 0U, 110U, 97U, 114U, 79U,
    63U, 33U, 23U, 0U, 0U, 44U, 39U, 39U, 0U, 0U, 39U, 0U, 44U, 72U, 80U, 41U,
    0U, 0U, 0U, 0U, 0U, 39U, 40U, 0U, 48U, 26U, 49U, 59U, 51U, 18U, 3U, 0U, 105U,
    73U, 0U, 0U, 0U, 0U, 60U, 60U, 0U, 0U, 0U, 13U, 0U, 0U, 23U, 54U, 0U, 0U,
    10U, 21U, 85U, 23U, 23U, 0U, 0U, 70U, 10U, 49U, 0U, 39U, 123U, 45U, 125U,
    109U, 86U, 77U, 117U, 58U, 47U, 47U, 9U, 83U, 96U, 0U, 0U, 0U, 0U, 14U, 39U,
    91U, 65U, 0U, 105U, 0U, 0U, 23U, 54U, 50U, 32U, 33U, 0U, 0U, 32U, 44U, 22U,
    30U, 37U, 39U, 36U, 39U, 0U, 0U, 72U, 0U, 80U, 10U, 26U, 0U, 10U, 18U, 40U,
    0U, 33U, 0U, 0U, 0U, 49U, 0U, 0U, 59U, 41U, 18U, 0U, 18U, 0U, 3U, 100U, 0U,
    0U, 49U, 0U, 32U, 60U, 33U, 0U, 0U, 0U, 0U, 25U, 26U, 2U, 0U, 18U, 0U, 41U,
    23U, 0U, 23U, 0U, 23U, 70U, 0U, 0U, 0U, 38U, 49U, 36U, 0U, 111U, 119U, 45U,
    3U, 125U, 106U, 109U, 88U, 27U, 86U, 10U, 25U, 76U, 41U, 58U, 23U, 10U, 0U,
    17U, 47U, 9U, 0U, 18U, 30U, 26U, 61U, 0U, 0U, 23U, 32U, 66U, 26U, 50U, 44U,
    105U, 51U, 0U, 0U, 0U, 18U, 23U, 0U, 0U, 0U, 0U, 0U, 32U, 0U, 0U, 0U, 0U,
    22U, 0U, 30U, 25U, 37U, 0U, 0U, 34U, 0U, 0U, 0U, 32U, 10U, 0U, 78U, 0U, 0U,
    0U, 0U, 0U, 0U, 0U, 0U, 40U, 40U, 33U, 10U, 21U, 0U, 59U, 10U, 0U, 18U, 18U,
    0U, 4U, 18U, 0U, 0U, 4U, 0U, 49U, 18U, 23U, 25U, 0U, 36U, 33U, 19U, 0U, 25U,
    0U, 0U, 0U, 2U, 0U, 0U, 10U, 0U, 23U, 0U, 0U, 0U, 23U, 23U, 0U, 21U, 0U, 38U,
    0U, 0U, 30U, 36U, 23U, 0U, 84U, 0U, 0U, 45U, 3U, 0U, 74U, 92U, 69U, 106U,
    107U, 0U, 0U, 0U, 9U, 27U, 25U, 0U, 0U, 0U, 25U, 0U, 0U, 76U, 41U, 0U, 21U,
    0U, 23U, 20U, 10U, 0U, 17U, 0U, 0U, 0U, 0U, 9U, 0U, 18U, 30U, 0U, 0U, 26U,
    61U, 0U, 0U, 0U, 21U, 0U, 0U, 18U, 0U, 26U, 32U, 36U, 0U, 10U, 0U, 105U, 18U,
    18U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 22U, 0U, 0U, 16U, 25U, 0U, 0U, 10U, 34U,
    18U, 0U, 0U, 0U, 36U, 23U, 39U, 0U, 23U, 0U, 0U, 18U, 0U, 0U, 0U, 21U, 57U,
    23U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 4U, 10U, 0U, 4U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
    25U, 12U, 36U, 0U, 0U, 0U, 19U, 3U, 0U, 0U, 2U, 0U, 0U, 0U, 0U, 0U, 0U, 10U,
    9U, 23U, 21U, 0U, 0U, 0U, 30U, 0U, 0U, 0U, 1U, 23U, 21U, 67U, 23U, 0U, 0U,
    0U, 38U, 74U, 92U, 90U, 69U, 69U, 106U, 106U, 91U, 0U, 0U, 9U, 0U, 0U, 0U,
    25U, 18U, 0U, 43U, 30U, 0U, 0U, 0U, 21U, 0U, 0U, 6U, 0U, 0U, 0U, 17U, 5U, 0U,
    0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 21U, 0U, 0U, 10U, 26U, 0U, 32U, 23U, 0U,
    0U, 0U, 0U, 10U, 0U, 0U, 0U, 0U, 0U, 22U, 0U, 16U, 11U, 0U, 0U, 0U, 0U, 34U,
    0U, 0U, 0U, 0U, 0U, 0U, 0U, 39U, 25U, 23U, 0U, 0U, 0U, 0U, 0U, 18U, 41U, 23U,
    0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 12U, 10U, 0U, 0U, 0U, 0U, 3U, 0U, 0U, 0U,
    0U, 0U, 9U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 1U, 0U, 0U, 0U, 0U, 21U, 46U, 67U,
    0U, 0U, 0U, 38U, 0U, 0U, 82U, 0U, 90U, 24U, 18U, 68U, 64U, 69U, 106U, 32U,
    0U, 81U, 0U, 23U, 0U, 0U, 0U, 0U, 0U, 0U, 4U, 43U, 30U, 0U, 0U, 0U, 6U, 0U,
    0U, 0U, 5U, 0U, 0U, 0U, 0U, 0U, 0U, 26U, 0U, 23U, 23U, 14U, 0U, 0U, 18U, 0U,
    0U, 0U, 11U, 0U, 32U, 34U, 18U, 39U, 0U, 25U, 0U, 23U, 0U, 0U, 0U, 41U, 0U,
    10U, 0U, 12U, 10U, 0U, 0U, 0U, 0U, 0U, 0U, 1U, 0U, 0U, 18U, 46U, 66U, 0U, 0U,
    0U, 23U, 25U, 33U, 56U, 0U, 24U, 0U, 18U, 0U, 0U, 55U, 64U, 0U, 0U, 89U,
    103U, 0U, 32U, 10U, 0U, 23U, 0U, 0U, 4U, 0U, 0U, 0U, 0U, 6U, 0U, 0U, 0U, 0U,
    0U, 0U, 0U, 23U, 0U, 14U, 9U, 18U, 0U, 0U, 0U, 32U, 26U, 34U, 0U, 0U, 0U, 0U,
    39U, 0U, 0U, 0U, 0U, 0U, 18U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 1U, 18U, 18U, 39U,
    0U, 44U, 0U, 0U, 15U, 0U, 25U, 33U, 10U, 0U, 56U, 0U, 24U, 0U, 0U, 55U, 29U,
    0U, 0U, 89U, 0U, 21U, 32U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
    0U, 0U, 8U, 9U, 0U, 0U, 32U, 0U, 0U, 10U, 0U, 10U, 0U, 38U, 0U, 0U, 0U, 0U,
    0U, 0U, 18U, 9U, 0U, 0U, 0U, 18U, 0U, 15U, 0U, 0U, 0U, 0U, 0U, 10U, 0U, 30U,
    0U, 0U, 0U, 0U, 29U, 0U, 89U, 89U, 0U, 21U, 0U, 0U, 8U, 0U, 0U, 0U, 0U, 0U,
    0U, 0U, 0U, 0U, 25U, 38U, 0U, 0U, 9U, 0U, 18U, 0U, 0U, 0U, 0U, 0U, 30U, 0U,
    0U, 29U, 56U, 0U, 0U, 18U, 0U, 0U, 8U, 3U, 0U, 25U, 32U, 10U, 0U, 0U, 0U, 0U,
    0U, 0U, 0U, 0U, 56U, 52U, 18U, 0U, 8U, 0U, 3U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
    10U, 55U, 52U, 0U, 0U, 0U, 0U, 8U, 0U, 0U, 0U, 0U, 0U, 55U, 0U, 0U, 0U, 0U,
    0U, 0U };

  int idx;
  double scores[29];
  double unusedU4[29];
  int k;
  boolean_T exitg1;
  m = 1.0;
  while ((m <= 1003.0) && (!(pruneList[(int)m - 1] <= 0)) && (!rtIsNaN(X[(int)
           obj_CutVar[(int)m - 1] - 1]))) {
    if (X[(int)obj_CutVar[(int)m - 1] - 1] < obj_CutPoint[(int)m - 1]) {
      m = obj_Children[((int)m - 1) << 1];
    } else {
      m = obj_Children[1 + (((int)m - 1) << 1)];
    }
  }

  for (idx = 0; idx < 29; idx++) {
    scores[idx] = obj_ClassProb[((int)m + 1003 * idx) - 1];
  }

  for (idx = 0; idx < 29; idx++) {
    unusedU4[idx] = 0.0;
    for (k = 0; k < 29; k++) {
      unusedU4[idx] += scores[k] * obj_Cost[k + 29 * idx];
    }
  }

  if (!rtIsNaN(unusedU4[0])) {
    idx = 1;
  } else {
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 30)) {
      if (!rtIsNaN(unusedU4[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (idx == 0) {
    idx = 1;
  } else {
    m = unusedU4[idx - 1];
    for (k = idx; k + 1 < 30; k++) {
      if (m > unusedU4[k]) {
        m = unusedU4[k];
        idx = k + 1;
      }
    }
  }

  return obj_ClassNames[idx - 1];
}

/*
 * File trailer for CompactClassificationTree.c
 *
 * [EOF]
 */
