#ifndef __EXPROJ_TYPEDEF_H__
#define __EXPROJ_TYPEDEF_H__

#include <atomic>

#include <dg_sparse/DG_Sparse.hpp>

#include "exproj_utils.h"

#include <stdexcept>


// #define CHECK_AVERAGE
#define ENABLE_CRIETERIA


#define MAX_ITER ( 1000 )
#define EPSILON  ( 1e-7 )

#define MAX_NNZ_PER_COL ( 10 )

#define D2R( x ) ( ( x ) * M_PI / 180. )
#define R2D( x ) ( ( x ) * 180. / M_PI )

#define IDX( x, y, n ) ( ( x ) + ( n ) * ( y ) )



namespace ExProj {
    enum n_ExProj
    {
        _n_stt = 7,
        _n_ctr = 4
    };

    struct admm_params {
        double sg1;
        double sg2;
        double sg3;
        double sg4;
        double sg5;

        double weight_f;
        double weight_p;
        double weight_v;
    };

    struct planet_data {
        double gravity;
    };

    struct lander_data {
        double m_wet;
        double m_dry;

        double T_max;
        double T_min;

        double mdot_pT;

        double max_theta;

        double z_lb[500];
        double z_ub[500];
        double s_lb[500];
        double s_ub[500];
    };

    struct trajec_data {
        double r_init[3];
        double v_init[3];
        double T_init[3];

        double r_des[3];
        double v_des[3];

        double max_gamma;
        double max_speed;
    };

    struct scale_param {
        double scale_p;
        double scale_v;
        double scale_a;
        double scale_s;
        double scale_kg;
        double scale_f;
    };
}

using namespace ExProj;

class SolverDG {
    private:
    /* ADMM parameters */
    const int max_N_TIM = 0;
    const int n_stt = _n_stt;
    const int n_ctr = _n_ctr;
    const int n_var_unit = n_stt + n_ctr;

    double t_pdg = -1;
    int    n_tim = -1;
    int    n_var = -1;
    int    n_eqc = -1;
    int    n_aug = -1;
    int    n_iq1 = -1;
    int    n_iq2 = -1;
    int    n_iq3 = -1;
    int    n_iq4 = -1;
    int    n_iq5 = -1;
    int    n_iqc = -1;

    int    n_iqc_0 = -1;
    int    n_iqc_1 = -1;
    int    n_iqc_2 = -1;
    int    n_iqc_3 = -1;
    int    n_iqc_4 = -1;

    struct admm_params aparam{};
    struct planet_data planet{};
    struct lander_data lander{};
    struct trajec_data trajec{};
    struct scale_param sparam{};

    bool scaled = false;

    public:
    /* LU decomposition */
    SparseLU<double>& LUdecomp;

    /* Sparse Matrices */
    SparseMatrix<double>& C_aug;
    SparseMatrix<double>& D_icn;
    SparseMatrix<double>& D_ict;

    /* Vectors */
    double* d_buf = nullptr;
    double* d_aug = nullptr;
    double* d_icn = nullptr;

    /* ADMM Variables */
    double* ynk = nullptr;
    double* unk = nullptr;
    double* znk = nullptr;
    double* snk = nullptr;

    double* ynp = nullptr;
    double* znp = nullptr;
    double* lnk = nullptr;
    double* lnp = nullptr;

    /* Buffers for C_aug triplets */
    int*    C_aug_rows = nullptr;
    int*    C_aug_cols = nullptr;
    double* C_aug_data = nullptr;

    /* Buffers for D_icx triplets */
    int*    D_icn_rows = nullptr;
    int*    D_icn_cols = nullptr;
    double* D_icn_data = nullptr;
    double* D_ict_data = nullptr;

    /* ***** Parameters ************* */
    // double T_min;
    // double T_max;

    // double delta_max;
    // double gamma_max;
    // double deltaTmax;

    // double gravity;
    // double mdot_pT;
    /* ****************************** */

    /* Indicates solver's status */
    std::atomic<bool> info{};
    std::atomic<int>  iter{};

    Eigen::Matrix<double,2,MAX_ITER> residual = Eigen::Matrix<double,2,MAX_ITER>::Zero();
    Eigen::VectorXd p_resid;
    Eigen::VectorXd d_resid;

    constexpr std::size_t get_n_tim() const { return n_tim; }
    constexpr std::size_t get_n_aug() const { return n_aug; }
    constexpr std::size_t get_n_var() const { return n_var; }
    constexpr std::size_t get_n_iqc() const { return n_iqc; }
    constexpr std::size_t get_n_var_unit() const { return n_var_unit; }

    inline SolverDG( SparseLU<double>& _LUdecomp, SparseMatrix<double>& _C_aug, SparseMatrix<double>& _D_icn, SparseMatrix<double>& _D_ict, const int _max_N_TIM );
    inline ~SolverDG();

    inline void update_param( admm_params& _aparam, lander_data& _lander, planet_data& _planet, trajec_data& _trajec, double _t_pdg, int _n_tim );

    void solver_setup( int verbose );
    void solver_solve( int verbose );

    inline void analyze_leq( const int save );
    inline void scaling_var( double scale_p, double scale_s, double scale_kg );

    inline void call_buffer( const int _n_aug, const int _n_iqc );
    inline void init_solver( void );
    inline void unscale_sol();

    inline void clear_all();
    inline void free_all();
};


SolverDG::SolverDG( SparseLU<double>& _LUdecomp, SparseMatrix<double>& _C_aug, SparseMatrix<double>& _D_icn, SparseMatrix<double>& _D_ict, const int _max_N_TIM )
: LUdecomp( _LUdecomp ), C_aug( _C_aug ), D_icn( _D_icn ), D_ict( _D_ict ), max_N_TIM( _max_N_TIM )
{
    const int n_var = n_var_unit * max_N_TIM;
    const int n_eqc = n_stt      * max_N_TIM;
    const int n_aug = n_var + n_eqc;

    const int n_iq1 = max_N_TIM * n_ctr;
    const int n_iq2 = max_N_TIM * 2 - 1;
    const int n_iq3 = max_N_TIM *   1  ;
    const int n_iq4 = max_N_TIM *   3  ;
    const int n_iq5 = max_N_TIM *   3  ;

    const int n_iqc_0 = 0;
    const int n_iqc_1 = n_iq1;
    const int n_iqc_2 = n_iq1 + n_iq2;
    const int n_iqc_3 = n_iq1 + n_iq2 + n_iq3;
    const int n_iqc_4 = n_iq1 + n_iq2 + n_iq3 + n_iq4;
    const int n_iqc_5 = n_iq1 + n_iq2 + n_iq3 + n_iq4 + n_iq5;
    const int n_iqc   = n_iq1 + n_iq2 + n_iq3 + n_iq4 + n_iq5;

    call_buffer( n_aug, n_iqc );

    D_icn.call_buffer();
    D_ict.call_buffer();

    p_resid = Eigen::VectorXd::Zero( n_iqc );
    d_resid = Eigen::VectorXd::Zero( n_iqc );
}


void SolverDG::call_buffer( const int _n_aug, const int _n_iqc )
{
    if ( C_aug_rows == nullptr ) { C_aug_rows = new int   [_n_aug*MAX_NNZ_PER_COL]; }
    if ( C_aug_cols == nullptr ) { C_aug_cols = new int   [_n_aug*MAX_NNZ_PER_COL]; }
    if ( C_aug_data == nullptr ) { C_aug_data = new double[_n_aug*MAX_NNZ_PER_COL]; }

    if ( D_icn_rows == nullptr ) { D_icn_rows = new int   [_n_iqc*MAX_NNZ_PER_COL]; }
    if ( D_icn_cols == nullptr ) { D_icn_cols = new int   [_n_iqc*MAX_NNZ_PER_COL]; }
    if ( D_icn_data == nullptr ) { D_icn_data = new double[_n_iqc*MAX_NNZ_PER_COL]; }
    if ( D_ict_data == nullptr ) { D_ict_data = new double[_n_iqc*MAX_NNZ_PER_COL]; }

    if ( d_buf == nullptr ) { d_buf = new double[_n_aug]; }
    if ( d_aug == nullptr ) { d_aug = new double[_n_aug]; }
    if ( d_icn == nullptr ) { d_icn = new double[_n_iqc]; }

    if ( ynk == nullptr ) { ynk = new double[_n_aug]; }
    if ( unk == nullptr ) { unk = new double[_n_iqc]; }
    if ( znk == nullptr ) { znk = new double[_n_iqc]; }
    if ( snk == nullptr ) { snk = new double[_n_iqc]; }

    if ( ynp == nullptr ) { ynp = new double[_n_aug]; }
    if ( znp == nullptr ) { znp = new double[_n_iqc]; }
    if ( lnk == nullptr ) { lnk = new double[_n_iqc]; }
    if ( lnp == nullptr ) { lnp = new double[_n_iqc]; }
}


void SolverDG::update_param( admm_params& _aparam, lander_data& _lander, planet_data& _planet, trajec_data& _trajec, double _t_pdg, int _n_tim )
{
    aparam = _aparam;
    lander = _lander;
    planet = _planet;
    trajec = _trajec;

    t_pdg = _t_pdg;
    n_tim = _n_tim;

    if ( n_tim > max_N_TIM )
    {
        throw std::runtime_error( "Future out of memory error, n_tim must smaller than max_N_TIM" );
    }
    else
    {
        n_var = n_var_unit * n_tim;
        n_eqc = n_stt      * n_tim;
        n_aug = n_var + n_eqc;

        n_iq1 = n_tim * n_ctr;
        n_iq2 = n_tim * 2 - 1;
        n_iq3 = n_tim *   1  ;
        n_iq4 = n_tim *   3  ;
        n_iq5 = n_tim *   3  ;

        n_iqc_0 = 0;
        n_iqc_1 = n_iq1;
        n_iqc_2 = n_iq1 + n_iq2;
        n_iqc_3 = n_iq1 + n_iq2 + n_iq3;
        n_iqc_4 = n_iq1 + n_iq2 + n_iq3 + n_iq4;
        n_iqc   = n_iq1 + n_iq2 + n_iq3 + n_iq4 + n_iq5;

        C_aug.reshape( n_aug, n_aug );
        D_icn.reshape( n_aug, n_aug );
        D_ict.reshape( n_aug, n_aug );
        LUdecomp.reshaped();
    }
}


SolverDG::~SolverDG()
{
    free_all();
}


void SolverDG::analyze_leq( const int save )
{
    solver_setup( 1 );

    Eigen::MatrixXd C_debug_d( n_aug, n_aug );
    C_aug.DenseView( C_debug_d );

    if ( save )
    {
        std::string filepath0 = "/home/jarvis/workspace/plot/datas/C_augs/C_aug" + std::to_string(int(n_aug)) + ".csv";

        save_csv( filepath0, C_debug_d );
    }

    Eigen::SparseMatrix<double> C_debug_s = C_debug_d.sparseView();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> EigenLU;

    EigenLU.analyzePattern( C_debug_s );

    const Eigen::VectorXi& Q = EigenLU.colsPermutation().indices();

    LUdecomp.ColumnOrder( Q.data() );
}


void SolverDG::scaling_var( double scale_p, double scale_s, double scale_kg )
{
    double scale_v = scale_p / scale_s;
    double scale_a = scale_v / scale_s;
    double scale_f = scale_a * scale_kg;

    sparam.scale_p  = scale_p;
    sparam.scale_v  = scale_v;
    sparam.scale_a  = scale_a;
    sparam.scale_s  = scale_s;
    sparam.scale_kg = scale_kg;
    sparam.scale_f  = scale_f;

    if ( !scaled )
    {
        t_pdg /= scale_s;

        lander.m_wet   /= scale_kg;
        lander.m_dry   /= scale_kg;
        lander.mdot_pT /= ( scale_s / scale_p );
        lander.T_max   /= scale_f;
        lander.T_min   /= scale_f;

        planet.gravity /= scale_a;

        trajec.r_init[0] /= scale_p; trajec.r_des[0] /= scale_p;
        trajec.r_init[1] /= scale_p; trajec.r_des[1] /= scale_p;
        trajec.r_init[2] /= scale_p; trajec.r_des[2] /= scale_p;

        trajec.v_init[0] /= scale_v; trajec.v_des[0] /= scale_v;
        trajec.v_init[1] /= scale_v; trajec.v_des[1] /= scale_v;
        trajec.v_init[2] /= scale_v; trajec.v_des[2] /= scale_v;

        trajec.max_speed /= scale_v;

        const double dt = t_pdg / (double)(n_tim);

        for ( int i = 0; i < n_tim; i++ )
        {
            lander.z_lb[i] = log( lander.m_wet - (i+1) * lander.mdot_pT * lander.T_max * dt );
            lander.z_ub[i] = log( lander.m_wet - (i+1) * lander.mdot_pT * lander.T_min * dt );

            lander.s_lb[i] = lander.T_min / exp( lander.z_ub[i] );
            lander.s_ub[i] = lander.T_max / exp( lander.z_lb[i] );
        }
    }

    scaled = true;
}


void SolverDG::init_solver( void )
{
    int idx;

    double dt    = t_pdg / (double)(n_tim);
    double coef  = 0.;
    double mag_T = lander.T_min;

    for ( int i = 0; i < n_tim; i++ )
    {
        idx  = n_var_unit * i;
        coef = (double)(i+1) / n_tim;

        ynk[idx+ 0] =  0.;
        ynk[idx+ 1] =  0.;
        ynk[idx+ 2] = -mag_T;
        ynk[idx+ 3] =  mag_T;

        ynk[idx+ 4] = ( trajec.r_des[0] - trajec.r_init[0] ) * coef + trajec.r_init[0];
        ynk[idx+ 5] = ( trajec.r_des[1] - trajec.r_init[1] ) * coef + trajec.r_init[1];
        ynk[idx+ 6] = ( trajec.r_des[2] - trajec.r_init[2] ) * coef + trajec.r_init[2];

        ynk[idx+ 7] = ( trajec.v_des[0] - trajec.v_init[0] ) * coef + trajec.v_init[0];
        ynk[idx+ 8] = ( trajec.v_des[1] - trajec.v_init[1] ) * coef + trajec.v_init[1];
        ynk[idx+ 9] = ( trajec.v_des[2] - trajec.v_init[2] ) * coef + trajec.v_init[2];

        ynk[idx+10] = log( lander.m_wet ) - abs( mag_T * lander.mdot_pT ) * coef * dt;
    }

    D_icn.gaxpy( znk, ynk );

    for ( int i = 0; i < n_iqc; i++ ) { unk[i] = znk[i]; }
}


void SolverDG::unscale_sol()
{
    int idx = 0;

    for ( int i = 0; i < ( n_tim + 1 ); i++ )
    {
        idx = n_var_unit * i;        

        ynk[idx+ 0] *= sparam.scale_a;
        ynk[idx+ 1] *= sparam.scale_a;
        ynk[idx+ 2] *= sparam.scale_a;
        ynk[idx+ 3] *= sparam.scale_a;

        ynk[idx+ 4] *= sparam.scale_p;
        ynk[idx+ 5] *= sparam.scale_p;
        ynk[idx+ 6] *= sparam.scale_p;

        ynk[idx+ 7] *= sparam.scale_v;
        ynk[idx+ 8] *= sparam.scale_v;
        ynk[idx+ 9] *= sparam.scale_v;
        
        ynk[idx+10] *= sparam.scale_kg;
    }
}


void SolverDG::clear_all()
{
    const int max_n_C_aug = n_aug * MAX_NNZ_PER_COL;
    const int max_n_D_icx = n_iqc * MAX_NNZ_PER_COL;

    for ( int i = 0; i < max_n_C_aug; i++ )
    {
        C_aug_rows[i] = 0;
        C_aug_cols[i] = 0;
        C_aug_data[i] = 0;
    }

    for ( int i = 0; i < max_n_D_icx; i++ )
    {
        D_icn_rows[i] = 0;
        D_icn_cols[i] = 0;
        D_icn_data[i] = 0;
        D_ict_data[i] = 0;
    }

    for ( int i = 0; i < n_aug; i++ )
    {
        d_buf[i] = 0;
        d_aug[i] = 0;
    }

    for ( int i = 0; i < n_aug; i++ )
    {
        ynk[i] = 0;
        ynp[i] = 0;
    }

    for ( int i = 0; i < n_iqc; i++ )
    {
        d_icn[i] = 0;

        unk[i] = 0;
        znk[i] = 0;
        snk[i] = 0;

        znp[i] = 0;
        lnk[i] = 0;
        lnp[i] = 0;
    }
}


void SolverDG::free_all()
{
    if ( C_aug_rows != nullptr ) delete[] C_aug_rows;
    if ( C_aug_cols != nullptr ) delete[] C_aug_cols;
    if ( C_aug_data != nullptr ) delete[] C_aug_data;

    if ( D_icn_rows != nullptr ) delete[] D_icn_rows;
    if ( D_icn_cols != nullptr ) delete[] D_icn_cols;
    if ( D_icn_data != nullptr ) delete[] D_icn_data;
    if ( D_ict_data != nullptr ) delete[] D_ict_data;

    if ( d_buf != nullptr ) delete[] d_buf;
    if ( d_aug != nullptr ) delete[] d_aug;
    if ( d_icn != nullptr ) delete[] d_icn;

    if ( ynk != nullptr ) delete[] ynk;
    if ( unk != nullptr ) delete[] unk;
    if ( znk != nullptr ) delete[] znk;
    if ( snk != nullptr ) delete[] snk;

    if ( ynp != nullptr ) delete[] ynp;
    if ( znp != nullptr ) delete[] znp;
    if ( lnk != nullptr ) delete[] lnk;
    if ( lnp != nullptr ) delete[] lnp;
}

#endif