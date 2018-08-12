#include "domain3D.h"

Domain_3::Domain_3(const Vec3d & l, const Vec3d & origin)
  : l_(l), half_l_(0.5 * l), origin_(origin), end_pnt_(origin + l) { 
}

ExtDomain_3::ExtDomain_3(const Vec3d & l, const Vec3d & origin,
                         const Vec_3<bool>& flag_ext)
  : Domain_3(l, origin), flag_ext_(flag_ext) {}


