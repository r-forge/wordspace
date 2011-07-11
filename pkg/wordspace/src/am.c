/*
 *  Compute different association measures from frequency signatures
 */

#include "function_defs.h"
#include "am.h"

double am_frequency(double f, double f1, double f2, double N, int sparse) {
  return f;
}

double am_simple_ll(double f, double f1, double f2, double N, int sparse) {
  double O = f, E = f1 * f2 / N;
  double ll = 2 * ( ((O > 0) ? O * log(O / E) : 0) - (O - E) );
  if (sparse)
    return (O > E) ? ll : 0;
  else
    return (O >= E) ? ll : -ll;
}

double am_t_score(double f, double f1, double f2, double N, int sparse) {
  double O = f, E = f1 * f2 / N;
  if (sparse) 
    return (O > E) ? (O - E) / sqrt(O) : 0;
  else
    return (O - E) / sqrt(O + 1); /* "discounted" t-score for O == 0 */
}

double am_z_score(double f, double f1, double f2, double N, int sparse) {
  double O = f, E = f1 * f2 / N;
  if (sparse)
    return (O > E) ? (O - E) / sqrt(E) : 0;
  else
    return (O - E) / sqrt(E); /* E == 0 should never happen */
}

double am_Dice(double f, double f1, double f2, double N, int sparse) {
  return 2 * f / (f1 + f2);
}

double am_MI(double f, double f1, double f2, double N, int sparse) {
  double O = f, E = f1 * f2 / N;
  if (sparse)
    return (O > E) ? log2(O / E) : 0;
  else
    return log2(O / E); /* not clear how to avoid the -Inf result here */
}

double transform(double x, int method) {
  switch (method) {
    case 0:       /* 0 = none */
      return x;
    case 1:       /* 1 = signed log */
      return sign(x) * log(fabs(x) + 1);
    case 2:       /* 2 = signed square root */
      return sign(x) * sqrt(fabs(x));
    case 3:       /* 3 = sigmoid (tanh) */
      return tanh(x);
    default:
      error("score transformation #%d is not defined -- internal error", method);
  }
}

int am_table_entries = 6;

am_func am_table[] = {
  &am_frequency,
  &am_simple_ll,
  &am_t_score,
  &am_z_score,
  &am_Dice,
  &am_MI
};

