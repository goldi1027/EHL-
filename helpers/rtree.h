#pragma once
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include "point.h"
#include "scenario.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace pl = polyanya;

namespace boost { namespace geometry { namespace traits {
  template<> struct tag<pl::Point> {
    typedef point_tag type;
  };

  template<> struct coordinate_type<pl::Point> {
    typedef double type;
  };

  template<> struct coordinate_system<pl::Point> {
    typedef cs::cartesian type;
  };

  template<> struct dimension<pl::Point> : boost::mpl::int_<2> {};

  template<> struct access<pl::Point, 0> {
    static double get(pl::Point const& p) {
      return p.x;
    }

    static void set(pl::Point& p, double const& value) {
      p.x = value;
    }
  };

  template<> struct access<pl::Point, 1> {
    static double get(pl::Point const& p) {
      return p.y;
    }

    static void set(pl::Point& p, double const& value) {
      p.y = value;
    }
  };

} } }

template <typename Rtree, typename Id>
size_t remove_ids_bulk(Rtree& rtree, Id const& id) {
    using V = typename Rtree::value_type;
    std::vector<V> v;
    std::copy_if(rtree.begin(), rtree.end(), back_inserter(v), [id](V const& v) { return v.second == id; });

    return rtree.remove(v.begin(), v.end());
}
