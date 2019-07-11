#include "material/material.h"
#include "material/bingham.h"

// Bingham 2D
static Register<mpm::Material<2>, mpm::Bingham<2>, unsigned, const Json&>
    bingham_2d("Bingham2D");

// Bingham 3D
static Register<mpm::Material<3>, mpm::Bingham<3>, unsigned, const Json&>
    bingham_3d("Bingham3D");
