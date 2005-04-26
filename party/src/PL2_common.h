
/* include R header files */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */

/* include private header files: this need to be restricted */

#include "Classes.h"
#include "Utils.h"
#include "mvt.h"
#include "LinearStatistic.h"
#include "TestStatistic.h"
#include "Distributions.h"
#include "Convenience.h"
#include "S3Classes.h"
#include "IndependenceTest.h"
#include "Splits.h"
#include "Node.h"
#include "Predict.h"
#include "SurrogateSplits.h"
