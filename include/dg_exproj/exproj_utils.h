#ifndef __EXPROJ_UTILS_H__
#define __EXPROJ_UTILS_H__

#include <iostream>
#include <string>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/Sparse>



#define EXPROJ_EPSILON ( 1e-5 )



double func( double z, const double z_ref, const double s_ref, const double param );
double grad( double z, const double z_ref, const double s_ref, const double param );

void proj_via_ExProj( double *x_opt, const double *x_ref, const double min_thrust, const double max_thrust );

Eigen::MatrixXd read_csv( const std::string& filePath );

void save_csv( const std::string& filePath, const Eigen::MatrixXd data );

#endif