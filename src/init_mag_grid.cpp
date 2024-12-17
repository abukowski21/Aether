// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// ----------------------------------------------------------------------
// Routine to convert p and q to r and theta. Can be solved iteratively,
// or with approach from (Swisdak, 2006), who solved it analytically:
//  https://arxiv.org/pdf/physics/0606044
//
// ----------------------------------------------------------------------
std::pair<precision_t, precision_t> qp_to_r_theta(precision_t q, precision_t p)
{
  // return quanties
  precision_t r, theta;
  // Intermediate quantities:
  precision_t term0, term1, term2, term3;

  term0 = 256.0 / 27.0 * pow(q, 2.0) * pow(p, 4.0);
  term1 = pow((1.0 + sqrt(1.0 + term0)), 2.0 / 3.0);
  term2 = pow(term0, 1.0 / 3.0);
  term3 = 0.5 * pow(((pow(term1, 2) + term1 * term2 + pow(term2, 2)) / term1), 3.0 / 2.0);

  r = p * (4.0 * term3) / (1.0 + term3) / (1.0 + sqrt(2.0 * term3 - 1.0));

  // now that r is determined we can solve for theta
  // theta = asin(sqrt(r/p));
  theta = asin(q * pow(r, 2.0));

  return {r, theta};
}

arma_vec baselat_spacing(precision_t extent,
                         precision_t origin,
                         precision_t upper_lim,
                         precision_t lower_lim,
                         int16_t nLats,
                         precision_t spacing_factor)
{
  std::string function = "Grid::baselat_spacing";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // intermediate latitude values
  precision_t lat_low, lat_high, lat_low0, lat_high0;
  // intermediate calculation values
  precision_t dlat, bb, aa, ang0, angq;

  // Now we can allocate the return array,
  arma_vec Lats(nLats);

  // Noting the special case of 1 root node & 1 processor...
  bool DO_FLIPBACK = false;
  if (extent > 0.5)
  {
    DO_FLIPBACK = true;
    nLats = nLats / 2;
    extent = 0.5;
  }

  // get the upper & lower latitude bounds for our division of the quadree
  if (origin < 0)
  {
    // negative origin: lat_high <=> lat_low,  & vice-versa.
    lat_low, lat_high = -upper_lim, -lower_lim;
    lat_low0 = lat_low;
    lat_low = lat_low - (lat_high - lat_low) * (origin / 0.5);
    lat_high = lat_low0 - (lat_high - lat_low0) * (extent / .5 + origin / 0.5);
  }
  else
  {
    lat_low0 = lower_lim;
    lat_low = lat_low + (lat_high - lat_low) * (origin / 0.5);
    lat_high = lat_low0 + (lat_high - lat_low0) * (extent / .5 + origin / 0.5);
  }

  // normalized spacing in latitude
  // NOTE: spacing factor != 1 will not work yet. but framework is here...
  bb = (lat_high - lat_low) / (pow(lat_high, spacing_factor) - pow(lat_low, spacing_factor));
  aa = lat_high - bb * pow(lat_high, spacing_factor);
  dlat = (lat_high - lat_low) / (nLats);
  report.print(4, "baselates laydown!");

  for (int64_t j = 0; j < nLats; j++)
  {
    ang0 = lat_low + (j + 1) * dlat;
    angq = aa + bb * pow(ang0, spacing_factor);
    Lats[j] = angq;
  }
  report.print(5, "baselates flipback!");

  if (DO_FLIPBACK)
    for (int64_t j = 0; j < nLats; j++)
    {
      Lats[j + nLats] = -1 * Lats[j];
    }
  report.print(4, "baselates flipbackdone!");

  report.exit(function);
  return Lats;
}

// // Gravity vectors in the dipole basis
// void calc_dipole_gravity(Planets planet){

// // rhat = -(2*cos/(del)) qhat + (sin/(del)) phat


// }


// === SPACING ALONG FIELD LINE === //
// Coordinates along the field line to begin modeling
// - In dipole (p,q) coordinates
// - North & south hemisphere base
  
void Grid::fill_field_lines(arma_vec baseLats, int64_t nAlts,
                            precision_t min_altRe, precision_t Gamma,
                            Planets planet)
{

  std::string function = "Grid::fill_field_lines";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nLats = baseLats.n_elem;

  precision_t q_S, q_N, delqp;
  arma_mat bAlts(nLats, nAlts), bLats(nLats, nAlts);

  // allocate & calculate some things outside of the main loop
  // - mostly just factors to make the code easier to read
  precision_t qp0, fb0, ft, delq, qp2, fa, fb, term0, term1, term2, term3, new_r;
  // exp_q_dist is the fraction of total q-distance to step for each pt along field line
  arma_vec exp_q_dist(nAlts);
  
  // temp holding of results from q,p -> r,theta conversion:
  std:: pair<precision_t, precision_t> r_theta;

  // Find L-Shell for each baseLat
  // using L=R/sin2(theta), where theta is from north pole
  arma_vec Lshells(nLats);
  for (int64_t iLat = 0; iLat < nLats; iLat++)
    Lshells(iLat) = (min_altRe) / pow(sin(cPI / 2 - baseLats(iLat)), 2.0);
  report.print(3, "lshells calculated!");
  
  // for (int64_t iLon = 0; iLon < nLons; iLon ++){
    for (int64_t iLat = 0; iLat < nLats; iLat ++){
      magP_scgc.col(iLat) = magP_scgc.col(iLat) + Lshells(iLat);
  }

  int64_t nZby2 = nAlts / 2;

  for (int64_t iAlt = 0; iAlt < nAlts; iAlt++)
    exp_q_dist(iAlt) = Gamma + (1 - Gamma) * exp(-pow(((iAlt - nZby2) / (nAlts / 10.0)), 2.0));
  report.print(3, "expQ");

  for (int iLat = 0; iLat < nLats; iLat++)
  {
    q_S = -cos(cPI / 2 + baseLats(iLat)) / pow(min_altRe, 2.0);
    q_N = -q_S;

    // calculate const. stride similar to sami2/3 (huba & joyce 2000) 
    // Note, this is not the:
    // ==  >>   sinh(gamma*qi)/sinh(gamma*q_S)  <<  ==
    // but a different formula where the spacing is easily controlled.
    // Doesn't have any lat/lon dependence so won't work for offset dipoles
    delqp = (q_N - q_S) / (nAlts + 1);
    delqp = min_altRe * delqp;
    for (int iAlt = 0; iAlt < nAlts; iAlt++)
    {
      qp0 = q_S + iAlt * (delqp);
      fb0 = (1 - exp_q_dist(iAlt)) / exp(-q_S / delqp - 1);
      ft = exp_q_dist(iAlt) - fb0 + fb0 * exp(-(qp0 - q_S) / delqp);
      delq = qp0 - q_S;

      // Q value at this point:
      qp2 = q_S + ft * delq;
      for (int64_t iLon=0; iLon < nLons; iLon ++) {
        magQ_scgc(iLon, iLat, iAlt) = qp2;
      }

      r_theta = qp_to_r_theta(qp2, Lshells(iLat));
      bAlts(iLat, iAlt) = r_theta.first;
      bLats(iLat, iAlt) = r_theta.second;
    }
  }

  report.print(3, "QP-rtheta done!");

  arma_vec rNorm1d(nAlts), lat1dAlong(nAlts);
  arma_cube r3d(nLons, nLats, nAlts);
  precision_t planetRadius;

  // rad_unit_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  // This is wrong, but get_radius doesnt support latitude at the time of writing
  planetRadius =  planet.get_radius(bLats(0)); 

  for (int64_t iLat = 0; iLat < nLats; iLat++)
  {
    for (int64_t iLon = 0; iLon < nLons; iLon++)
    {
      // Not currently used. Dipole isn't offset. Leaving just in case.
      // Lon = magPhi_scgc(iLon, iLat, 1);

      for (int64_t iAlt = 0; iAlt < nAlts; iAlt++)
      {
        lat1dAlong(iAlt) = bLats(iLat, iAlt);
        rNorm1d(iAlt) = bAlts(iLat, iAlt);
      }
      // Indexing is weird, but consistent. Might be helpful to plot out...
      r3d.tube(iLon,  nLats - iLat - 1) = rNorm1d * planetRadius;
      magLat_scgc.tube(iLon, nLats - iLat - 1) = lat1dAlong;
    }
  }

  magAlt_scgc = r3d;

  report.exit(function);
  return;
}

////////////////////////////////////////////
// convert cell coordinates to geographic //
////////////////////////////////////////////
std::vector <arma_cube> mag_to_geo(arma_cube magLon, arma_cube magLat, arma_cube magAlt,
Planets planet){
  std::string function = "Grid::init_dipole_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::vector<arma_cube> llr, xyz_mag, xyz_geo, xyzRot1, xyzRot2;
  llr.push_back(magLon);
  llr.push_back(magLat);
  llr.push_back(magAlt);
  xyz_mag = transform_llr_to_xyz_3d(llr);

  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  std::vector<precision_t> dipole_center = planet.get_dipole_center();

  // Reverse our dipole rotations:
  xyzRot1 = rotate_around_y_3d(xyz_mag, magnetic_pole_tilt);
  xyzRot2 = rotate_around_z_3d(xyzRot1, magnetic_pole_rotation);

  // offset dipole (not yet implemented):
  // xyz_geo[0] = xyzRot2[0] + dipole_center[0];
  // xyz_geo[1] = xyzRot2[1] + dipole_center[1];
  // xyz_geo[2] = xyzRot2[2] + dipole_center[2];

  // transform back to lon, lat, radius:
  llr = transform_xyz_to_llr_3d(xyzRot2);
  return llr;
}

// -----------------------------------------------------------------------
// Convert XyzDipole to XyzGeo
//  
// -----------------------------------------------------------------------

void Grid::convert_dipole_geo_xyz(Planets planet, precision_t XyzDipole[3], precision_t XyzGeo[3]) {
  
  std::string function = "Grid::init_dipole_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t XyzRemoveShift[3];
  precision_t XyzRemoveTilt[3];
  precision_t XyzRemoveRot[3];

  // get planetary parameters
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t radius = planet.get_radius(0.0);

  
  // get the dipole shift, but normalize it to equatorial radius 
  precision_t dipole_center[3];
  std::vector<precision_t> temp_dipole_center = planet.get_dipole_center();
  if ((temp_dipole_center[0] != 0) or (temp_dipole_center[1] != 0) or (temp_dipole_center[2] != 0)){
    report.print(0, "Dipole center != 0, but that is not supported yet. Setting to 0!");
    temp_dipole_center = {0,0,0};

  }

  transform_float_vector_to_array(temp_dipole_center, dipole_center);

  dipole_center[0]=dipole_center[0]/radius;
  dipole_center[1]=dipole_center[1]/radius;
  dipole_center[2]=dipole_center[2]/radius;

  // Remove Tilt
  transform_rot_y(XyzDipole, magnetic_pole_tilt, XyzRemoveTilt);

  // Remove Rot
  transform_rot_z(XyzRemoveTilt, magnetic_pole_rotation, XyzRemoveRot);

  // Remove Shift
  vector_add(XyzRemoveRot, dipole_center, XyzGeo);  
}

// ----------------------------------------------------------------------
// Initialize the dipole grid.
// - inputs (min_apex, min_alt, LatStretch, FieldLineStretch, max_lat_dipole)
//   are read from input files. And the numbers of each coordinate.
// - nLats must be even!!
// ----------------------------------------------------------------------
bool Grid::init_dipole_grid(Quadtree quadtree_ion, Planets planet)
{

  using namespace std;
  bool DidWork = true;


  string function = "Grid::init_dipole_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // turn the switch on!
  IsGeoGrid = false;
  IsMagGrid = true;
  IsCubeSphereGrid=false;

  report.print(0, "Creating inter-node connections Grid");
  if (!Is0D & !Is1Dz) create_sphere_connection(quadtree_ion);

  report.print(0, "Creating Dipole Grid");

  report.print(3, "Getting mgrid_inputs inputs in dipole grid");

  Inputs::grid_input_struct grid_input = input.get_grid_inputs("ionGrid");

  // Number of ghost cells:
  int64_t nGCs = get_nGCs();

  // Get inputs:
  precision_t min_alt = grid_input.alt_min * cKMtoM;
  precision_t LatStretch = grid_input.LatStretch;
  precision_t Gamma = grid_input.FieldLineStretch;
  precision_t min_apex = grid_input.min_apex * cKMtoM;
  precision_t max_lat = grid_input.max_blat;

  // Normalize inputs to planet radius... (update when earth is oblate)
  precision_t planetRadius = planet.get_radius(0.0);
  // Altitude to begin modeling, normalized to planet radius
  precision_t min_alt_re = (min_alt + planetRadius) / planetRadius;
  precision_t min_apex_re = (min_apex + planetRadius) / planetRadius;

  if (LatStretch != 1){
    report.error("LatStretch values =/= 1 are not yet supported!");
    DidWork=false;}

  if (nAlts % 2 != 0){
    report.error("nAlts must be even!");
    DidWork=false;}

  if (min_alt >= min_apex){
    report.error("min_apex must be more than min_alt");
    DidWork=false;}

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree_ion.get_vect("LL"); // origin
  arma_vec size_right_norm = quadtree_ion.get_vect("SR"); // lon_lims
  arma_vec size_up_norm = quadtree_ion.get_vect("SU");    //[1] = lat_lims
  report.print(3, "longitudes");

  precision_t dlon = size_right_norm(0) * cPI / (nLons - 2 * nGCs);
  precision_t lon0 = lower_left_norm(0) * cPI;
  arma_vec lon1d(nLons);
  // if we are not doing anything in the lon direction, then set dlon to
  // something reasonable:
  if (!HasXdim)
    dlon = 1.0 * cDtoR;

  // Dimension iterators
  int64_t iLon, iLat, iAlt;

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      magLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
      // geoLon_Left.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
  }

  report.print(3, "longitudes point two");

  // Latitudes:

  // min_lat calculated from min_apex
  precision_t min_lat = acos(sqrt(1 / min_apex_re));

  // latitude of field line base:
  // todo: needs support for variable stretching. it's like, halfway there.
  arma_vec baseLats = baselat_spacing(size_up_norm[1], lower_left_norm[1],
                                      max_lat, min_lat, nLats, 1.0);

  report.print(3, "baselats done!");

  // latitude & altitude of points on field lines (2D)
  fill_field_lines(baseLats, nAlts, min_apex_re, Gamma, planet);
  
  
  report.print(3, "Done generating symmetric latitude & altitude spacing in dipole.");


  std::vector <arma_cube> llr = mag_to_geo(magLon_scgc, magLat_scgc, magAlt_scgc, planet);

  geoLon_scgc = llr[0];
  geoLat_scgc = llr[1];
  geoAlt_scgc = llr[2] - planetRadius;
  report.print(4, "Done dipole -> geographic transformations for the dipole grid.");

  // Calculate the radius, of planet
  fill_grid_radius(planet);

  // Figure out what direction is radial:
  rad_unit_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  gravity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iV = 0; iV < 3; iV++)
  {
    rad_unit_vcgc[iV].zeros();
    gravity_vcgc[iV].zeros();
  }

  arma_cube br = 2 * sin(abs(magLat_scgc));
  arma_cube bt = cos(magLat_scgc);
  arma_cube bm = sqrt(br % br + bt % bt);
  // Latitudinal direction of radial:
  arma_cube s = sign(magLat_scgc);
  // s.elem(find(s == 0)).ones();

  rad_unit_vcgc[1] = bt / bm % s;
  rad_unit_vcgc[2] = -br / bm;

  precision_t mu = planet.get_mu();
  gravity_vcgc[1] = mu * rad_unit_vcgc[1] % radius2i_scgc;
  gravity_vcgc[2] = mu * rad_unit_vcgc[2] % radius2i_scgc;
  gravity_potential_scgc.set_size(nX, nY, nAlts);
  gravity_potential_scgc.zeros();
  gravity_mag_scgc = sqrt(
      gravity_vcgc[0] % gravity_vcgc[0] +
      gravity_vcgc[1] % gravity_vcgc[1] +
      gravity_vcgc[2] % gravity_vcgc[2]);

  report.print(4, "Done gravity calculations for the dipole grid.");

  calc_dipole_grid_spacing(planet);
  report.print(4, "Done altitude spacing for the dipole grid.");

  // Calculate magnetic field and magnetic coordinates:
  fill_grid_bfield(planet);
  report.print(4, "Done filling dipole grid with b-field!");


  // put back into altitude. we've been carrying around radius:
  magAlt_scgc = magAlt_scgc - planetRadius;

  report.exit(function);
  return DidWork;
}
