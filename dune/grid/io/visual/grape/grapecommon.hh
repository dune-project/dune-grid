// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRAPECOMMON_HH
#define DUNE_GRAPECOMMON_HH

#include "grapewrapper.hh"

#if HAVE_GRAPE

/* max number for vertices and faces in Grape */
enum { MAX_EL_DOF  = 8 };
enum { MAX_EL_FACE = 6 };

/* global variables for maxlevel use */
static BUTTON * maxlevelButton=0;

/* on click set min and max value of function to colorbar */
static BUTTON * minMaxColorbar=0;

/* global variables for iterator choice */
static COMBOBUTTON  * iteratorButton = 0;
static int defaultIteratorValue = 0 ;

/* global variables for partition type choice */
static COMBOBUTTON * partitionTypeButton = 0;

static TIMESCENE * globalTsc = 0;

void setupLeafButton(MANAGER *mgr, void *sc, int yesTimeScene);
void removeLeafButton(MANAGER *mgr, void *sc);
void setDefaultIteratorValue(int val);

#endif // end HAVE_GRAPE

/* info about data on one mesh */
typedef struct datainfo DATAINFO;
struct datainfo
{
  const char * name;
  const char * base_name;
  DATAINFO *next;

  int dimVal; /* length of vector (dimVal = 1 --> scalar, otherwise vector  */
  int * comp; /* number of each component */
};

/* info about one mesh */
typedef struct info INFO;
struct info
{
  int fix_mesh; /* if no dynamic grid 1 : else 0 */
  const char  *name;
  DATAINFO * datinf;
  void  *tsc;
};
#endif // end DUNE_GRAPECOMMON_HH
