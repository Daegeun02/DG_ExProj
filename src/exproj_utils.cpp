#include "dg_exproj/exproj_utils.h"

#include <fstream>
#include <vector>
#include <iomanip>

#include <math.h>



double func( double z, const double z_ref, const double s_ref, const double param )
{
    double output = 0;
    double exp_z  = exp( z );

    output += exp_z * ( z - z_ref );
    output -= param * param / exp_z;
    output += param * s_ref;

    return output;
}


double grad( double z, const double z_ref, const double s_ref, const double param )
{
    double output = 0;
    double exp_z  = exp( z );

    output += exp_z * ( z - z_ref + 1.0 );
    output += param * param / exp_z;

    return output;

}


void proj_via_ExProj( double *x_opt, const double *x_ref, const double min_thrust, const double max_thrust )
{
    double z_opt = 0.00;
    double _func = 1.e8;
    double _grad = 1.e8;

    uint32_t max_iter = 100;
    uint32_t iter = 0;

    if ( x_ref[1] < ( min_thrust / exp( x_ref[0] ) ) )
    {
        while ( fabs( _func ) > (double)(EXPROJ_EPSILON) )
        {
            _func = func( z_opt, x_ref[0], x_ref[1], min_thrust );
            _grad = grad( z_opt, x_ref[0], x_ref[1], min_thrust );

            z_opt -= _func / _grad;

            iter++;

            if ( iter == max_iter ) break;
        }

        x_opt[0] = z_opt;
        x_opt[1] = min_thrust / exp( z_opt );
    }
    else if ( x_ref[1] > ( max_thrust / exp( x_ref[0] ) ) )
    {
        while ( fabs( _func ) > (double)(EXPROJ_EPSILON) )
        {
            _func = func( z_opt, x_ref[0], x_ref[1], max_thrust );
            _grad = grad( z_opt, x_ref[0], x_ref[1], max_thrust );

            z_opt -= _func / _grad;

            iter++;

            if ( iter == max_iter ) break;
        }

        x_opt[0] = z_opt;
        x_opt[1] = max_thrust / exp( z_opt );
    }
    else
    {
        x_opt[0] = x_ref[0];
        x_opt[1] = x_ref[1];
    }
}


// Function to read CSV file into an Eigen matrix
Eigen::MatrixXd read_csv( const std::string& filePath ) 
{
    std::ifstream file(filePath);
    std::vector<std::vector<double>> data;
    std::string line;

    // Read data line by line
    while ( std::getline( file, line ) ) 
    {
        std::stringstream lineStream( line );
        std::string cell;
        std::vector<double> rowData;

        // Split each line by commas
        while ( std::getline( lineStream, cell, ',' ) ) 
        {
            rowData.push_back( std::stod( cell ) );
        }
        data.push_back( rowData );
    }

    // Determine the number of rows and columns
    int numRows = data.size();
    int numCols = ( (numRows > 0) ? data[0].size() : 0 );

    // Create an Eigen matrix
    Eigen::MatrixXd matrix( numRows, numCols );

    // Populate the Eigen matrix with data
    for ( int i = 0; i < numRows; ++i ) 
    {
        for ( int j = 0; j < numCols; ++j )
        {
            matrix(i, j) = data[i][j];
        }
    }
    return matrix;
}


void save_csv( const std::string& filePath, const Eigen::MatrixXd data )
{
    std::ofstream file( filePath );

    file << std::fixed << std::setprecision(10);

    for ( int j = 0; j < data.cols(); j++ )
    {
        file << j << ", ";
    }

    file << 0 << std::endl;

    for ( int i = 0; i < data.rows(); i++ )
    {
        for ( int j = 0; j < data.cols(); j++ )
        {
            file << data(i,j) << ", ";
        }

        file << 0 << std::endl;
    }

    file.close();
}