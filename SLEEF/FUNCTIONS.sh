#!/bin/bash
set -x
./FUNCTION.sh sin10 sin35 cos10 cos35 tan10 tan35 \
              log10 log35 log210 log235 log1p10 log1010 exp10 exp210 exp1010 expm110 \
              sqrt05 sqrt35 cbrt10 cbrt35 \
              asin10 asin35 acos10 acos35 atan10 atan35 \
              sinh10 sinh35 cosh10 cosh35 tanh10 tanh35 \
              asinh10 acosh10 atanh10 \
              erf10 erfc15 tgamma10 lgamma10
./FUNCTION_1_2.sh atan210 atan235 pow10 hypot05 hypot35
./FUNCTION_2_1.sh sincos10 sincos35
