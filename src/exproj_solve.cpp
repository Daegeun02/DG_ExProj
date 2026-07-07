#include "dg_exproj/exproj_typedef.h"



void SolverDG::solver_solve( int verbose )
{
    for ( int i = 0; i < MAX_ITER; i++ )
    {
        residual(0,i) = 1e8;
        residual(1,i) = 1e8;
    }

    const double dt = t_pdg / (double)(n_tim);

    // init_solver();

    double buf[4];
    double value;
    double temp0;
    double temp1;
    double temp2;

    const double T_min = lander.T_min;
    const double T_max = lander.T_max;

    const double gimbal = cos( lander.max_theta );
    const double glides = tan( trajec.max_gamma );
    const double maxvel = trajec.max_speed;

    const double* z_lb = lander.z_lb;
    const double* z_ub = lander.z_ub;
    const double* s_lb = lander.s_lb;
    const double* s_ub = lander.s_ub;

    double* ptr_p_resid = p_resid.data();
    double* ptr_d_resid = d_resid.data();

    uint32_t idx1 = 0;
    uint32_t idx2 = 0;

    info = false;
    iter = -1;

    for ( int _iter = 0; _iter < MAX_ITER; _iter++ )
    {
        for ( int i = 0; i < n_iqc; i++ ) { znp[i] = znk[i]; };
        
        /* ****************************************************************** */
        /* ************** znk & lnk update BEGIN **************************** */
        /* ****************************************************************** */
        for ( int i = 0; i < n_tim; i++ )
        {
            /* ************************************************************** */
            /* ********** D1w, LossLess Contraint BEGIN ********************* */
            /* ************************************************************** */
            idx1 = i * n_var_unit;
            idx2 = n_iqc_0 + n_ctr * i;

            buf[0] = unk[idx2+0] + lnk[idx2+0];
            buf[1] = unk[idx2+1] + lnk[idx2+1];
            buf[2] = unk[idx2+2] + lnk[idx2+2];
            buf[3] = unk[idx2+3] + lnk[idx2+3];

            value = sqrt(
                buf[0]*buf[0] + buf[1]*buf[1] + buf[2]*buf[2]
            );

            if ( value != buf[3] )
            {
                if ( buf[3] > 0 )
                {
                    buf[3] += value;
                    buf[3] *= 0.5;

                    value  /= buf[3];

                    buf[0] /= value;
                    buf[1] /= value;
                    buf[2] /= value;
                }
                else
                {
                    buf[0] = 0.;
                    buf[1] = 0.;
                    buf[2] = 0.;
                    buf[3] = 0.;
                }
            }

            znk[idx2+0] = buf[0];
            znk[idx2+1] = buf[1];
            znk[idx2+2] = buf[2];
            znk[idx2+3] = buf[3];
            /* ************************************************************** */
            /* ********** D1w, LossLess Contraint  END  ********************* */
            /* ************************************************************** */

            /* ************************************************************** */
            /* ********** D2w, ExProj Thrust Contraint BEGIN **************** */
            /* ************************************************************** */
            idx1 = i * n_var_unit;
            idx2 = n_iqc_1 + 2 * i;

            if ( i != 0 )
            {
                buf[0] = unk[idx2-1] + lnk[idx2-1];
                buf[1] = unk[idx2+0] + lnk[idx2+0];

                if ( buf[1] > s_ub[i] ) buf[1] = s_ub[i];
                if ( buf[1] < s_lb[i] ) buf[1] = s_lb[i];

                proj_via_ExProj( &(buf[2]), &(buf[0]), T_min, T_max );

                znk[idx2-1] = buf[2];
                znk[idx2+0] = buf[3];
            }
            else
            {
                buf[0] = log( lander.m_wet );
                buf[1] = unk[idx2+0] + lnk[idx2+0];

                if ( buf[1] > s_ub[i] ) buf[1] = s_ub[i];
                if ( buf[1] < s_lb[i] ) buf[1] = s_lb[i];

                proj_via_ExProj( &(buf[2]), &(buf[0]), T_min, T_max );

                znk[idx2+0] = buf[3];
            }
            /* ************************************************************** */
            /* ********** D2w, ExProj Thrust Contraint  END  **************** */
            /* ************************************************************** */

            /* ************************************************************** */
            /* ********** D3w, Gimbal Limit Contraint BEGIN ***************** */
            /* ************************************************************** */
            idx1 = i * n_var_unit;
            idx2 = n_iqc_2 + i;

            buf[2] = unk[idx2+0] + lnk[idx2+0];

            if ( buf[2] > 0 ) buf[2] = 0.;

            znk[idx2+0] = buf[2];
            /* ************************************************************** */
            /* ********** D3w, Gimbal Limit Contraint  END  ***************** */
            /* ************************************************************** */

            /* ************************************************************** */
            /* ********** D4w, Glide Slope Contraint BEGIN ****************** */
            /* ************************************************************** */
            idx1 = i * n_var_unit;
            idx2 = n_iqc_3 + 3 * i;

            buf[0] = unk[idx2+0] + lnk[idx2+0];
            buf[1] = unk[idx2+1] + lnk[idx2+1];
            buf[2] = unk[idx2+2] + lnk[idx2+2];

            value = sqrt( buf[0]*buf[0] + buf[1]*buf[1] );

            if ( buf[2] < ( value / glides ) )
            {
                if ( buf[2] > 0 )
                {
                    temp2 = glides * glides;
                    temp0 = ( temp2 * value + glides * buf[2] ) / ( 1 + temp2 );
                    temp1 = temp0 / value;

                    buf[0] *= temp1;
                    buf[1] *= temp1;
                    buf[2]  = temp0 / glides;
                }
                else
                {
                    buf[0] = 0.;
                    buf[1] = 0.;
                    buf[2] = 0.;
                }
            }

            znk[idx2+0] = buf[0];
            znk[idx2+1] = buf[1];
            znk[idx2+2] = buf[2];
            /* ************************************************************** */
            /* ********** D4w, Glide Slope Contraint  END  ****************** */
            /* ************************************************************** */

            /* ************************************************************** */
            /* ********** D5w, Max Velocity Contraint BEGIN ***************** */
            /* ************************************************************** */
            idx1 = i * n_var_unit;
            idx2 = n_iqc_4 + 3 * i;

            buf[0] = unk[idx2+0] + lnk[idx2+0];
            buf[1] = unk[idx2+1] + lnk[idx2+1];
            buf[2] = unk[idx2+2] + lnk[idx2+2];

            value = sqrt(
                buf[0]*buf[0] + buf[1]*buf[1] + buf[2]*buf[2]
            );

            if ( value > maxvel )
            {
                temp0 = maxvel / value;

                buf[0] *= temp0;
                buf[1] *= temp0;
                buf[2] *= temp0;
            }

            znk[idx2+0] = buf[0];
            znk[idx2+1] = buf[1];
            znk[idx2+2] = buf[2];
            /* ************************************************************** */
            /* ********** D5w, Max Velocity Contraint  END  ***************** */
            /* ************************************************************** */
        }
        /* ****************************************************************** */
        /* ************** znk & lnk update  END  **************************** */
        /* ****************************************************************** */

        for ( int i = 0; i < n_iqc; i++ )
        {
            snk[i] = znk[i] - lnk[i];
            unk[i] = 0;
        }

        /* ****************************************************************** */
        /* ************** ynk update BEGIN ********************************** */
        /* ****************************************************************** */
        for ( int i = 0; i < n_aug; i++ ) { d_buf[i] = d_aug[i]; };

        D_ict.gaxpy( d_buf, snk );

        LUdecomp.Solve( d_buf );

        for ( int i = 0; i < n_aug; i++ ) { ynk[i] = d_buf[i]; };

        D_icn.gaxpy( unk, ynk );
        /* ****************************************************************** */
        /* ************** ynk update  END  ********************************** */
        /* ****************************************************************** */

        for ( int i = 0; i < n_iqc; i++ )
        {
            ptr_d_resid[i] = znk[i] - znp[i];
            ptr_p_resid[i] = unk[i] - znk[i];

            lnk[i] += ptr_p_resid[i];
        }

        /* ****************************************************************** */
        /* ************** Calculate Residual BEGIN ************************** */
        /* ****************************************************************** */
        residual(0,_iter) = d_resid.squaredNorm(); // / n_tim;
        residual(1,_iter) = p_resid.squaredNorm(); // / n_tim;

        // if ( ( residual(0,iter) < EPSILON ) && ( residual(1,iter) < EPSILON ) )
#ifdef ENABLE_CRIETERIA
        if ( residual(1,_iter) < EPSILON )
        {
            info =  true;
            iter = _iter;

            if ( verbose )
            {
                std::cout << std::fixed << std::setprecision(15);
                std::cout << "Solver solved at: " << iter+1 << " iteration" << std::endl;
                std::cout << "Dual Residual is  : " << residual(0,_iter) << std::endl;
                std::cout << "Primal Residual is: " << residual(1,_iter) << std::endl;
            }

            // std::cout << std::endl;
            // std::cout << d_resid.block<N_IQ1,1>(N_IQC_0,0).squaredNorm() << std::endl;
            // std::cout << d_resid.block<N_IQ2,1>(N_IQC_1,0).squaredNorm() << std::endl;
            // std::cout << d_resid.block<N_IQ3,1>(N_IQC_2,0).squaredNorm() << std::endl;
            // std::cout << p_resid.block<N_IQ4,1>(N_IQC_3,0).squaredNorm() << std::endl;
            // std::cout << p_resid.block<N_IQ5,1>(N_IQC_4,0).squaredNorm() << std::endl;
            // std::cout << p_resid.block<N_IQ6,1>(N_IQC_5,0).squaredNorm() << std::endl;
            // std::cout << std::endl;
            break;
        }
#endif
        /* ****************************************************************** */
        /* ************** Calculate Residual  END  ************************** */
        /* ****************************************************************** */
    }

    if ( info != true )
    {
        std::cout << std::fixed << std::setprecision(20);
        std::cout << "Solver failed to find solution, " << iter << std::endl;
        std::cout << "Dual Residual is  : " << residual(0,MAX_ITER-1) << std::endl;
        std::cout << "Primal Residual is: " << residual(1,MAX_ITER-1) << std::endl;

        if ( verbose )
        {
            std::cout << d_resid.block( n_iqc_0, 0, n_iq1, 1 ).squaredNorm() << std::endl;
            std::cout << d_resid.block( n_iqc_1, 0, n_iq2, 1 ).squaredNorm() << std::endl;
            std::cout << d_resid.block( n_iqc_2, 0, n_iq3, 1 ).squaredNorm() << std::endl;
            // std::cout << d_resid.block( N_IQC_3, 0, N_IQ4, 1 ).squaredNorm() << std::endl;
            // std::cout << d_resid.block( N_IQC_4, 0, N_IQ5, 1 ).squaredNorm() << std::endl;
            // std::cout << d_resid.block( N_IQC_5, 0, N_IQ6, 1 ).squaredNorm() << std::endl;
            std::cout << std::endl;
        }
    }
}