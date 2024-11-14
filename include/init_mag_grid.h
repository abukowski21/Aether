// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md


#ifndef INCLUDE_INIT_MAG_GRID_H_
#define INCLUDE_INIT_MAG_GRID_H_

#include "aether.h"
#include "planets.h"
#include "grid.h"

bool init_dipole_grid(Grid &mGrid, Planets planet);

// Analytic solution to get from q,p dipole coords to r,theta
// q coordinate along b-field line
// p l-shell 
// return (r,theta)
std::pair<precision_t, precision_t> qp_to_r_theta(precision_t q, precision_t p);

// get the latitude spacing given the quadtree start&size, and the latitude limits
// extent quadtree up
// origin quadtree origin
// upper_lim upper latitude limit (input)
// lower_lim lower latitude limit (from min_apex)
// nLats number of latitudes (nY)
// spacing_factor (not supported yet), so always 1.0.
arma_vec baselat_spacing(precision_t extent,
                        precision_t origin,
                        precision_t upper_lim,
                        precision_t lower_lim,
                        int16_t nLats,
                        precision_t spacing_factor);

// Take limits & specs of field line, fill it with points.
// std:: tuple <arma_mat, arma_mat, arma_vec> Grid::fill_field_lines (arma_vec lShells, 
//                                                 int64_t nAlts,
//                                                 int64_t nLats,
//                                                 precision_t Gamma);


// convert mag to geographic
std::vector <arma_cube> mag_to_geo(arma_cube magLon, 
                                   arma_cube magLat,
                                   arma_cube magAlt,
                                   Planets planet);


#endif  // INCLUDE_INIT_GEO_GRID_H_

