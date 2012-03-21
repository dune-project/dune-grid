// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __GRAPE_COMMON_CC__
#define __GRAPE_COMMON_CC__

#define MINIMUM(a,b) (((a) > (b)) ? (b) : (a))
#define MAXIMUM(a,b) (((a) < (b)) ? (b) : (a))

#include "grapecommon.hh"

enum { numberOfPartitionTypes = 6 };
static const char *partitionNames[numberOfPartitionTypes] = {
  "Interior_Partition",
  "InteriorBorder_Partition",
  "Overlap_Partition",
  "OverlapFront_Partition",
  "All_Partition",
  "Ghost_Partition"
};

enum { numberOfIterators = 4 };
static const char *iteratorNames[numberOfIterators] = {
  "LeafIterator",
  "LevelIterator",
  "Macro + HierarchicIterator",
  "GridPart Iterator",
};

static int actualPartitionType = 0;

static BUTTON *button_set_current_data_item(int pnr)
{
  BUTTON *lbutton = 0;

  lbutton = (BUTTON *)START_METHOD(G_INSTANCE);
  assert( lbutton );
  actualPartitionType = pnr;

  END_METHOD(lbutton);
}

/* set default iterator value */
inline void setDefaultIteratorValue(int val)
{
  assert( val >= 0 );
  assert( val < 4 );
  defaultIteratorValue = val;
}

/* add Button which can switch between LevelIteration and LeafIteration */
inline void setupLeafButton(MANAGER *mgr, void *sc, int yesTimeScene)
{
  assert(!iteratorButton);
  assert(!maxlevelButton);


  // partition types
  {
    int num = numberOfPartitionTypes;
    CYCLE_LABEL *clabel;

    clabel = (CYCLE_LABEL *) malloc(sizeof(CYCLE_LABEL) * (num + 1));

    for (int i=0; i<num; i++)
    {
      clabel[i].value = i;
      clabel[i].label = (char *) partitionNames[i];
    }
    // last of list is 0,and NULL
    clabel[num].value = 0;
    clabel[num].label = NULL;

    if(! (GRAPE_CALL(Button,"find-method") ("set-current-data-item")) )
    {
      GRAPE_CALL(Button,"add-method") ("set-current-data-item",
                                       button_set_current_data_item);
    }

    partitionTypeButton = (COMBOBUTTON *)GRAPE_CALL(ComboButton,"new-instance")
                            ("set-current-data-item",NULL,"",clabel);
    GRAPE_CALL(partitionTypeButton, "set-fill-mode") (MENU_FILL_BOTTOM);

    GRAPE_CALL(partitionTypeButton,"set-instance") (partitionTypeButton);
    GRAPE_CALL(partitionTypeButton,"set-label") (clabel[0].label);
    GRAPE_CALL(partitionTypeButton,"set-pref-size") (12.0,1.0);

    GRAPE_CALL(mgr,"add-inter") (partitionTypeButton);
  }

  // iterator types
  {
    int num = numberOfIterators;
    CYCLE_LABEL *clabel;

    clabel = (CYCLE_LABEL *) malloc(sizeof(CYCLE_LABEL) * (num + 1));

    for (int i=0; i<num; ++i)
    {
      clabel[i].value = i;
      clabel[i].label = (char *) iteratorNames[i];
    }
    // last of list is 0,and NULL
    clabel[num].value = 0;
    clabel[num].label = NULL;

    //GRAPE_CALL(Button,"add-method")("set-current-data-item",
    //                           button_set_current_data_item);

    iteratorButton = (COMBOBUTTON *)GRAPE_CALL(ComboButton,"new-instance")
                       ("set-current-data-item",NULL,"",clabel);
    GRAPE_CALL(iteratorButton, "set-fill-mode") (MENU_FILL_BOTTOM);

    GRAPE_CALL(iteratorButton,"set-instance") (iteratorButton);
    GRAPE_CALL(iteratorButton,"set-label") (clabel[0].label);
    GRAPE_CALL(iteratorButton,"set-pref-size") (12.0,1.0);

    GRAPE_CALL(iteratorButton,"set-value") (defaultIteratorValue);

    GRAPE_CALL(mgr,"add-inter") (iteratorButton);
  }

  maxlevelButton = (BUTTON *)
                   new_item (Button,
                             I_Label, "use only maxlevel",
                             I_Instance, sc,
                             I_Method, "maxlevel-on-off",
                             I_Size, 12., 1.,
                             I_FillMode, MENU_FILL_BOTTOM,
                             I_End);

  minMaxColorbar = (BUTTON *)
                   new_item (Button,
                             I_Label, "set min/max to Colorbar",
                             I_Instance, sc,
                             I_Method, "set-min-max-values",
                             I_Size, 12., 1.,
                             I_FillMode, MENU_FILL_BOTTOM,
                             I_End);

  minMaxColorbar = (BUTTON *)
                   new_item (Button,
                             I_Label, "set min/max to Colorbar",
                             I_Instance, sc,
                             I_Method, "set-min-max-values",
                             I_Size, 12., 1.,
                             I_FillMode, MENU_FILL_BOTTOM,
                             I_End);

  GRAPE_CALL(mgr,"add-inter") (maxlevelButton);
  GRAPE_CALL(mgr,"add-inter") (minMaxColorbar);

  GRAPE_CALL(maxlevelButton,"set-state") (PRESSED);
  GRAPE_CALL(minMaxColorbar,"set-state") (UNPRESSED);
  maxlevelButton->on_off = OFF;
}

/* add Button which can switch between LevelIteration and LeafIteration */
inline void removeLeafButton(MANAGER *mgr, void *sc)
{
  GRAPE_CALL(mgr,"remove-inter") (partitionTypeButton);
  GRAPE_CALL(partitionTypeButton,"delete") ();
  partitionTypeButton = 0;
  GRAPE_CALL(mgr,"remove-inter") (iteratorButton);
  GRAPE_CALL(iteratorButton,"delete") ();
  iteratorButton = 0;
  GRAPE_CALL(mgr,"remove-inter") (maxlevelButton);
  GRAPE_CALL(maxlevelButton,"delete") ();
  maxlevelButton = 0;
  GRAPE_CALL(mgr,"remove-inter") (minMaxColorbar);
  GRAPE_CALL(minMaxColorbar,"delete") ();
  minMaxColorbar = 0;
}

inline void timeSceneInit(INFO *info, const int n_info,
                          const int procs)
{
  // set global max partition number for partition display
  Dune :: __MaxPartition = procs;

  printf("Warning: set proc to 1\n");

  for (int n = 0; n < MAXIMUM(1, n_info); n++)
  {
    printf("n = %d, make TimeScene \n",n);
    info[n].tsc = (TIMESCENE *)GRAPE_CALL(TimeScene,"new-instance") (info[n].name);
    ((TIMESCENE *)info[n].tsc)->sync = 1;
    if (n > 0)
      ((TIMESCENE *) info[n-1].tsc)->next_scene = (SCENE *) info[n].tsc;
  }

  for (int n = 0; n < MAXIMUM(1, n_info); n++)
  {
    {
      TIMESCENE * tsc = (TIMESCENE *) info[n].tsc;
      // > 0 because tsc for proc 0 already exists
      for (int p = procs-1; p > 0; p--)
      {
        TIMESCENE * newSc = NULL;
        printf("add timescene for proc %d \n",p);
        assert(tsc);
        char * newName = (char *) malloc(strlen(info[n].name) + 5 * sizeof(char));
        assert(newName);
        sprintf(newName,"%s_%d",info[n].name,p);
        newSc = (TIMESCENE *)GRAPE_CALL(TimeScene,"new-instance") (newName);
        assert(newSc);

        newSc->sync = 1;
        newSc->next_scene = tsc->next_scene;
        tsc->next_scene = (SCENE *) newSc;
      }
    }

  }

  {
    TIMESCENE * newSc = NULL;
    printf("add timescene for proc combined object\n");
    char * newName = (char *) malloc(12 * sizeof(char));
    assert(newName);
    sprintf(newName,"combo obj");
    newSc = (TIMESCENE *)GRAPE_CALL(TimeScene,"new-instance") (newName);
    assert(newSc);

    newSc->sync = 1;
    newSc->next_scene = NULL;
    globalTsc = newSc;
  }
  return;
}

#if 0
/* add scene with combined object at the end of scene tree */
inline SCENE * combine_scenes_send ()
{
  SCENE* sc = (SCENE*) START_METHOD (G_INSTANCE);
  TIMESCENE * newSc = NULL;
  MANAGER * mgr = NULL;
  GRAPE_ALERT( sc, (char *) "combine-scenes: No hmesh!", END_METHOD(NULL));

  newSc = (TIMESCENE *) GRAPE_CALL(TimeScene,"new-instance") ("combined Scene");
  assert(newSc);

  while(sc)
  {
    if(sc->next_scene)
      sc = sc->next_scene;
    else
      break;
  }

  mgr = (MANAGER *)GRAPE_CALL(Manager,"get-stdmgr") ();
  assert(mgr);
  assert(globalTsc);

  sc->next_scene = (SCENE *) globalTsc;
  GRAPE_CALL(mgr, "goto-instance") (globalTsc);

  END_METHOD(sc);
}

/* call handle for a bunch of timescenes */
inline void displayTimeScene ( INFO *info, int procs )
{
  TIMESCENE *tsc = (TIMESCENE *) info[0].tsc;
  if(tsc)
  {
    MANAGER       *mgr;
#ifdef GRID_MODE
    GRAPHICDEVICE *grdev;

    grdev = (GRAPHICDEVICE *)GRAPE_CALL(GraphicDevice,"get-stddev") ();
    grdev->clear();
    if (grdev && (grdev->grid_patch != G_GRID))
    {
      GRAPE_CALL(grdev,"grid-patch") (G_GRID);
    }
#endif
    GrapeInterface_two_two::grape_add_remove_methods();
    GrapeInterface_two_three::grape_add_remove_methods();
    GrapeInterface_three_three::grape_add_remove_methods();
    GrapeInterface_three_three::initPartitionDisp( procs-1 );

    {
      /* add combine methods to send method of Scene and TimeScene */
      GRAPE_CALL(Scene,"add-method") ("combine-scenes-send",combine_scenes_send);
      GRAPE_CALL(TimeScene,"add-method") ("combine-scenes-send",combine_scenes_send);
    }

    mgr = (MANAGER *)GRAPE_CALL(Manager,"get-stdmgr") ();

    if(!maxlevelButton && !iteratorButton) setupLeafButton(mgr,tsc,1);

    GRAPE_CALL(mgr,"handle") (tsc);
  }
}
#endif


#undef MINIMUM
#undef MAXIMUM

#endif // #ifndef __GRAPE_COMMON_CC__
