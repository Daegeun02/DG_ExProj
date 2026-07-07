#include "dg_exproj/exproj_typedef.h"



void SolverDG::solver_setup( int verbose )
{
    const double dt = t_pdg / (double)(n_tim);

    uint32_t idx1;
    uint32_t idx2;
    uint32_t idx3;

    uint32_t C_aug_nnzs = 0;
    uint32_t D_icn_nnzs = 0;
    uint32_t D_ict_nnzs = 0;

    const double sg1 = aparam.sg1;
    const double sg2 = aparam.sg2;
    const double sg3 = aparam.sg3;
    const double sg4 = aparam.sg4;
    const double sg5 = aparam.sg5;

    const double gimbal = cos( lander.max_theta );

    const double value01 = (-0.5) * dt * dt;
    const double value02 = ( 1.0) * dt * lander.mdot_pT;
    const double value03 = ( 0.5) * planet.gravity * dt * dt;
    const double value04 = ( 1.0) * planet.gravity * dt;

    clear_all();

    for ( int i = 0; i < n_tim; i++ )
    {
        idx1 = i * n_stt + n_var;
        idx2 = i * n_var_unit;

        /* ****************************************************************** */
        /* ********** D1w, LossLess Constraint BEGIN ************************ */
        /* ****************************************************************** */
        idx3 = n_iqc_0 + n_ctr * i;

        D_icn_rows[D_icn_nnzs] = idx3+0; D_icn_cols[D_icn_nnzs] = idx2+0; D_icn_data[D_icn_nnzs++] = 1.0;
        D_icn_rows[D_icn_nnzs] = idx3+1; D_icn_cols[D_icn_nnzs] = idx2+1; D_icn_data[D_icn_nnzs++] = 1.0;
        D_icn_rows[D_icn_nnzs] = idx3+2; D_icn_cols[D_icn_nnzs] = idx2+2; D_icn_data[D_icn_nnzs++] = 1.0;
        D_icn_rows[D_icn_nnzs] = idx3+3; D_icn_cols[D_icn_nnzs] = idx2+3; D_icn_data[D_icn_nnzs++] = 1.0;

        D_ict_data[D_ict_nnzs++] = sg1;
        D_ict_data[D_ict_nnzs++] = sg1;
        D_ict_data[D_ict_nnzs++] = sg1;
        D_ict_data[D_ict_nnzs++] = sg1;
        /* ****************************************************************** */
        /* ********** D1w, LossLess Constraint  END  ************************ */
        /* ****************************************************************** */

        /* ****************************************************************** */
        /* ********** D2w, ExProj Thrust Constraint BEGIN ******************* */
        /* ****************************************************************** */
        idx3 = n_iqc_1 + 2 * i;

        D_icn_rows[D_icn_nnzs  ] = idx3+0; D_icn_cols[D_icn_nnzs] = idx2+3; D_icn_data[D_icn_nnzs++] = 1.0;
        D_ict_data[D_ict_nnzs++] = sg2;

        if ( i != 0 )
        {
            D_icn_rows[D_icn_nnzs  ] = idx3-1; D_icn_cols[D_icn_nnzs] = idx2-1; D_icn_data[D_icn_nnzs++] = 1.0;
            D_ict_data[D_ict_nnzs++] = sg2;
        }
        /* ****************************************************************** */
        /* ********** D2w, ExProj Thrust Constraint  END  ******************* */
        /* ****************************************************************** */

        /* ****************************************************************** */
        /* ********** D3w, Gimbal Limit Constraint BEGIN ******************** */
        /* ****************************************************************** */
        idx3 = n_iqc_2 + i;

        D_icn_rows[D_icn_nnzs] = idx3+0; D_icn_cols[D_icn_nnzs] = idx2+2; D_icn_data[D_icn_nnzs++] = 1.0;
        D_icn_rows[D_icn_nnzs] = idx3+0; D_icn_cols[D_icn_nnzs] = idx2+3; D_icn_data[D_icn_nnzs++] = gimbal;

        D_ict_data[D_ict_nnzs++] = sg3;
        D_ict_data[D_ict_nnzs++] = sg3 * gimbal;
        /* ****************************************************************** */
        /* ********** D3w, Gimbal Limit Constraint  END  ******************** */
        /* ****************************************************************** */

        /* ****************************************************************** */
        /* ********** D4w, Glide Slope Constraint BEGIN ********************* */
        /* ****************************************************************** */
        idx3 = n_iqc_3 + 3 * i;

        D_icn_rows[D_icn_nnzs] = idx3+0; D_icn_cols[D_icn_nnzs] = idx2+4; D_icn_data[D_icn_nnzs++] =  1.0;
        D_icn_rows[D_icn_nnzs] = idx3+1; D_icn_cols[D_icn_nnzs] = idx2+5; D_icn_data[D_icn_nnzs++] =  1.0;
        D_icn_rows[D_icn_nnzs] = idx3+2; D_icn_cols[D_icn_nnzs] = idx2+6; D_icn_data[D_icn_nnzs++] = -1.0;

        D_ict_data[D_ict_nnzs++] =  sg4;
        D_ict_data[D_ict_nnzs++] =  sg4;
        D_ict_data[D_ict_nnzs++] = -sg4;
        /* ****************************************************************** */
        /* ********** D4w, Glide Slope Constraint  END  ********************* */
        /* ****************************************************************** */

        /* ****************************************************************** */
        /* ********** D5w, Max Velocity Constraint BEGIN ******************** */
        /* ****************************************************************** */
        idx3 = n_iqc_4 + 3 * i;

        D_icn_rows[D_icn_nnzs] = idx3+0; D_icn_cols[D_icn_nnzs] = idx2+7; D_icn_data[D_icn_nnzs++] = 1.0;
        D_icn_rows[D_icn_nnzs] = idx3+1; D_icn_cols[D_icn_nnzs] = idx2+8; D_icn_data[D_icn_nnzs++] = 1.0;
        D_icn_rows[D_icn_nnzs] = idx3+2; D_icn_cols[D_icn_nnzs] = idx2+9; D_icn_data[D_icn_nnzs++] = 1.0;

        D_ict_data[D_ict_nnzs++] = sg5;
        D_ict_data[D_ict_nnzs++] = sg5;
        D_ict_data[D_ict_nnzs++] = sg5;
        /* ****************************************************************** */
        /* ********** D5w, Max Velocity Constraint  END  ******************** */
        /* ****************************************************************** */
    }

    D_icn.from_triplet( D_icn_rows, D_icn_cols, D_icn_data, D_icn_nnzs, 1 );
    D_ict.from_triplet( D_icn_cols, D_icn_rows, D_ict_data, D_ict_nnzs, 1 );

    D_ict.gabpc( D_icn, C_aug );

    /* ********************************************************************** */
    /* ********** Add Terminal state cost BEGIN ***************************** */
    /* ********************************************************************** */
    idx2 = ( n_tim - 1 ) * n_var_unit + 3;

    C_aug.x_ptr[C_aug.p_ptr[idx2+0]] += aparam.weight_f;
    C_aug.x_ptr[C_aug.p_ptr[idx2+1]] += aparam.weight_p;
    C_aug.x_ptr[C_aug.p_ptr[idx2+2]] += aparam.weight_p;
    C_aug.x_ptr[C_aug.p_ptr[idx2+3]] += aparam.weight_p;
    C_aug.x_ptr[C_aug.p_ptr[idx2+4]] += aparam.weight_v;
    C_aug.x_ptr[C_aug.p_ptr[idx2+5]] += aparam.weight_v;
    C_aug.x_ptr[C_aug.p_ptr[idx2+6]] += aparam.weight_v;
    /* ********************************************************************** */
    /* ********** Add Terminal state cost  END  ***************************** */
    /* ********************************************************************** */

    C_aug_nnzs = C_aug.to_triplet( C_aug_rows, C_aug_cols, C_aug_data );

    for ( int i = 0; i < n_tim; i++ )
    {
        idx1 = i * n_stt + n_var;
        idx2 = i * n_var_unit;

        /* ****************************************************************** */
        /* ********** Insert Input Matrix BEGIN ***************************** */
        /* ****************************************************************** */
        C_aug_rows[C_aug_nnzs] = idx1+0; C_aug_cols[C_aug_nnzs] = idx2+0; C_aug_data[C_aug_nnzs++] =  value01;
        C_aug_rows[C_aug_nnzs] = idx1+1; C_aug_cols[C_aug_nnzs] = idx2+1; C_aug_data[C_aug_nnzs++] =  value01;
        C_aug_rows[C_aug_nnzs] = idx1+2; C_aug_cols[C_aug_nnzs] = idx2+2; C_aug_data[C_aug_nnzs++] =  value01;
        C_aug_rows[C_aug_nnzs] = idx1+3; C_aug_cols[C_aug_nnzs] = idx2+0; C_aug_data[C_aug_nnzs++] = -dt;
        C_aug_rows[C_aug_nnzs] = idx1+4; C_aug_cols[C_aug_nnzs] = idx2+1; C_aug_data[C_aug_nnzs++] = -dt;
        C_aug_rows[C_aug_nnzs] = idx1+5; C_aug_cols[C_aug_nnzs] = idx2+2; C_aug_data[C_aug_nnzs++] = -dt;
        C_aug_rows[C_aug_nnzs] = idx1+6; C_aug_cols[C_aug_nnzs] = idx2+3; C_aug_data[C_aug_nnzs++] =  value02;

        C_aug_rows[C_aug_nnzs] = idx2+0; C_aug_cols[C_aug_nnzs] = idx1+0; C_aug_data[C_aug_nnzs++] =  value01;
        C_aug_rows[C_aug_nnzs] = idx2+1; C_aug_cols[C_aug_nnzs] = idx1+1; C_aug_data[C_aug_nnzs++] =  value01;
        C_aug_rows[C_aug_nnzs] = idx2+2; C_aug_cols[C_aug_nnzs] = idx1+2; C_aug_data[C_aug_nnzs++] =  value01;
        C_aug_rows[C_aug_nnzs] = idx2+0; C_aug_cols[C_aug_nnzs] = idx1+3; C_aug_data[C_aug_nnzs++] = -dt;
        C_aug_rows[C_aug_nnzs] = idx2+1; C_aug_cols[C_aug_nnzs] = idx1+4; C_aug_data[C_aug_nnzs++] = -dt;
        C_aug_rows[C_aug_nnzs] = idx2+2; C_aug_cols[C_aug_nnzs] = idx1+5; C_aug_data[C_aug_nnzs++] = -dt;
        C_aug_rows[C_aug_nnzs] = idx2+3; C_aug_cols[C_aug_nnzs] = idx1+6; C_aug_data[C_aug_nnzs++] =  value02;
        /* ****************************************************************** */
        /* ********** Insert Input Matrix  END  ***************************** */
        /* ****************************************************************** */

        /* ****************************************************************** */
        /* ********** Insert Identity Matrix BEGIN ************************** */
        /* ****************************************************************** */
        C_aug_rows[C_aug_nnzs] = idx1+ 0; C_aug_cols[C_aug_nnzs] = idx2+ 4; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx1+ 1; C_aug_cols[C_aug_nnzs] = idx2+ 5; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx1+ 2; C_aug_cols[C_aug_nnzs] = idx2+ 6; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx1+ 3; C_aug_cols[C_aug_nnzs] = idx2+ 7; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx1+ 4; C_aug_cols[C_aug_nnzs] = idx2+ 8; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx1+ 5; C_aug_cols[C_aug_nnzs] = idx2+ 9; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx1+ 6; C_aug_cols[C_aug_nnzs] = idx2+10; C_aug_data[C_aug_nnzs++] = 1.0;

        C_aug_rows[C_aug_nnzs] = idx2+ 4; C_aug_cols[C_aug_nnzs] = idx1+ 0; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx2+ 5; C_aug_cols[C_aug_nnzs] = idx1+ 1; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx2+ 6; C_aug_cols[C_aug_nnzs] = idx1+ 2; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx2+ 7; C_aug_cols[C_aug_nnzs] = idx1+ 3; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx2+ 8; C_aug_cols[C_aug_nnzs] = idx1+ 4; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx2+ 9; C_aug_cols[C_aug_nnzs] = idx1+ 5; C_aug_data[C_aug_nnzs++] = 1.0;
        C_aug_rows[C_aug_nnzs] = idx2+10; C_aug_cols[C_aug_nnzs] = idx1+ 6; C_aug_data[C_aug_nnzs++] = 1.0;
        /* ****************************************************************** */
        /* ********** Insert Identity Matrix  END  ************************** */
        /* ****************************************************************** */

        if ( i != 0 )
        {
            /* ************************************************************** */
            /* ********** Insert State Transition Matrix BEGIN ************** */
            /* ************************************************************** */
            C_aug_rows[C_aug_nnzs] = idx1+0; C_aug_cols[C_aug_nnzs] = idx2-7; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx1+1; C_aug_cols[C_aug_nnzs] = idx2-6; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx1+2; C_aug_cols[C_aug_nnzs] = idx2-5; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx1+3; C_aug_cols[C_aug_nnzs] = idx2-4; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx1+4; C_aug_cols[C_aug_nnzs] = idx2-3; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx1+5; C_aug_cols[C_aug_nnzs] = idx2-2; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx1+6; C_aug_cols[C_aug_nnzs] = idx2-1; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx1+0; C_aug_cols[C_aug_nnzs] = idx2-4; C_aug_data[C_aug_nnzs++] = -dt;
            C_aug_rows[C_aug_nnzs] = idx1+1; C_aug_cols[C_aug_nnzs] = idx2-3; C_aug_data[C_aug_nnzs++] = -dt;
            C_aug_rows[C_aug_nnzs] = idx1+2; C_aug_cols[C_aug_nnzs] = idx2-2; C_aug_data[C_aug_nnzs++] = -dt;

            C_aug_rows[C_aug_nnzs] = idx2-7; C_aug_cols[C_aug_nnzs] = idx1+0; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx2-6; C_aug_cols[C_aug_nnzs] = idx1+1; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx2-5; C_aug_cols[C_aug_nnzs] = idx1+2; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx2-4; C_aug_cols[C_aug_nnzs] = idx1+3; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx2-3; C_aug_cols[C_aug_nnzs] = idx1+4; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx2-2; C_aug_cols[C_aug_nnzs] = idx1+5; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx2-1; C_aug_cols[C_aug_nnzs] = idx1+6; C_aug_data[C_aug_nnzs++] = -1.0;
            C_aug_rows[C_aug_nnzs] = idx2-4; C_aug_cols[C_aug_nnzs] = idx1+0; C_aug_data[C_aug_nnzs++] = -dt;
            C_aug_rows[C_aug_nnzs] = idx2-3; C_aug_cols[C_aug_nnzs] = idx1+1; C_aug_data[C_aug_nnzs++] = -dt;
            C_aug_rows[C_aug_nnzs] = idx2-2; C_aug_cols[C_aug_nnzs] = idx1+2; C_aug_data[C_aug_nnzs++] = -dt;
            /* ************************************************************** */
            /* ********** Insert State Transition Matrix  END  ************** */
            /* ************************************************************** */

            /* ************************************************************** */
            /* ********** Insert Dynamic Constraint Vector BEGIN1 *********** */
            /* ************************************************************** */
            d_aug[idx1+2] = value03;
            d_aug[idx1+5] = value04;
            /* ************************************************************** */
            /* ********** Insert Dynamic Constraint Vector  END1  *********** */
            /* ************************************************************** */
        }
        else
        {
            /* ************************************************************** */
            /* ********** Insert Dynamic Constraint Vector BEGIN2 *********** */
            /* ************************************************************** */
            d_aug[idx1+0] = trajec.r_init[0] + trajec.v_init[0] * dt;
            d_aug[idx1+1] = trajec.r_init[1] + trajec.v_init[1] * dt;
            d_aug[idx1+2] = trajec.r_init[2] + trajec.v_init[2] * dt + value03;
            d_aug[idx1+3] = trajec.v_init[0];
            d_aug[idx1+4] = trajec.v_init[1];
            d_aug[idx1+5] = trajec.v_init[2] + value04;
            d_aug[idx1+6] = log( lander.m_wet );
            /* ************************************************************** */
            /* ********** Insert Dynamic Constraint Vector  END2  *********** */
            /* ************************************************************** */

            /* ************************************************************** */
            /* ********** Insert Cost Function Vector BEGIN ***************** */
            /* ************************************************************** */
            for ( int i = 0; i < n_tim; i++ )
            {
                idx1 = n_var_unit * i;

                d_aug[idx1+3] -= 1.0;
            }

            d_aug[n_var-8] += aparam.weight_f * planet.gravity;
            d_aug[n_var-7] += aparam.weight_p * trajec.r_des[0];
            d_aug[n_var-6] += aparam.weight_p * trajec.r_des[1];
            d_aug[n_var-5] += aparam.weight_p * trajec.r_des[2];
            d_aug[n_var-4] += aparam.weight_v * trajec.v_des[0];
            d_aug[n_var-3] += aparam.weight_v * trajec.v_des[1];
            d_aug[n_var-2] += aparam.weight_v * trajec.v_des[2];
            /* ************************************************************** */
            /* ********** Insert Cost Function Vector  END  ***************** */
            /* ************************************************************** */
        }
    }

    C_aug.from_triplet( C_aug_rows, C_aug_cols, C_aug_data, C_aug_nnzs, 1 );

    LUdecomp.Factorize();

    if ( verbose )
    {
        std::cout << "nnzs" << std::endl;
        std::cout << C_aug.nnzs() << std::endl;
        std::cout << LUdecomp.L.nnzs() << std::endl;
        std::cout << LUdecomp.U.nnzs() << std::endl;
    }
}