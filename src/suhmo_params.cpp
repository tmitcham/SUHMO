#include <suhmo_params.H>
#include "ParmParse.H"

#include "NamespaceHeader.H"

void suhmo_params::setDefaults()
{
    m_rho_i = 910.0;
    m_rho_w = 1000.0;
    m_gravity = 9.8;
    // constants of problem
    m_G = 0.0;
    m_L = 0.0;   
    m_H = 0.0;  
    m_nu = 0.0;
    m_ct = 0.0;
    m_cw = 0.0;
    m_omega = 0.0;
    m_basal_friction = true;
    // BC
    m_br = 0.0 ;
    m_lr = 0.0;
    m_A = 0.0;
    m_cutOffbr = 0.0;
    m_maxOffbr = 1000.0;
    m_DiffFactor = 0.0;
    // Initialization 
    m_slope = 0.0;
    m_gapInit = 0.0;
    m_ReInit = 0.0;
    // Moulin
    m_n_moulins = 0;
    m_distributed_input  = 0.0;
    m_time_varying_input = false;
    m_runoff             = 0.0;
    m_deltaT             = 0.0;

    m_ramp = false;
    m_ramp_up = 1.0;  
    m_duration_max = 1.0;
    m_relax =   0.1;
    m_floor_min = 0.0;
    m_floor_max = 1.0;

    m_coupled_to_bisicles = false;
    m_bisicles_input_file = "";
    m_output_N_file = "";
    m_ice_thickness_name = "thickness";
    m_bed_elevation_name = "Z_base";
    m_velocity_x_name = "xVel";
    m_velocity_y_name = "yVel";
}

void suhmo_params::readInputs()
{
    pout() << "suhmo_params::readInputs()" << endl;

    // Problem dependent variables
    ParmParse ppParams("suhmo");

    // Cst of problem
    ppParams.get("rho_i", m_rho_i);
    ppParams.get("rho_w", m_rho_w);
    ppParams.get("g", m_gravity);
    
    ppParams.get("GeoFlux", m_G);
    ppParams.get("LatHeat", m_L);
    ppParams.get("IceHeight", m_H);
    ppParams.get("WaterViscosity", m_nu);
    ppParams.get("ct", m_ct);
    ppParams.get("cw", m_cw);
    ppParams.get("turbulentParam", m_omega);
    ppParams.get("basalFriction", m_basal_friction);
    m_ub.resize(SpaceDim,0.0);
    ppParams.getarr("SlidingVelocity", m_ub, 0, SpaceDim);
    ppParams.get("br", m_br);
    ppParams.get("cutOffbr", m_cutOffbr);
    ppParams.get("maxOffbr", m_maxOffbr);
    ppParams.get("diffFactor", m_DiffFactor);
    ppParams.get("lr", m_lr);
    ppParams.get("A", m_A);
    ppParams.get("slope", m_slope);
    ppParams.get("GapInit", m_gapInit);
    ppParams.get("ReInit", m_ReInit);
    // MOULINS
    ppParams.get("time_varying_input", m_time_varying_input);
    ppParams.get("ramp", m_ramp);
    ppParams.get("distributed_input", m_distributed_input);
    ppParams.get("n_moulins", m_n_moulins);
    if (m_n_moulins > 0 ) {
        m_moulin_position.resize(m_n_moulins*SpaceDim,0.0);
        ppParams.getarr("moulin_position", m_moulin_position, 0, m_n_moulins*SpaceDim);
        m_moulin_flux.resize(m_n_moulins,0.0);
        ppParams.getarr("moulin_flux", m_moulin_flux, 0, m_n_moulins);
        m_sigma.resize(m_n_moulins,0.0);
        ppParams.getarr("moulin_sigma", m_sigma, 0, m_n_moulins);
        if (m_time_varying_input) {
            ppParams.get("Ra", m_runoff);
        }
        if (m_ramp) {
            ppParams.get("ramp_up", m_ramp_up);  
            ppParams.get("duration_max", m_duration_max);
            ppParams.get("relax", m_relax);
            ppParams.get("floor_min", m_floor_min);
            ppParams.get("floor_max", m_floor_max);
        }
    } else if (m_n_moulins < 0) {
        if (m_time_varying_input) {
            ppParams.get("deltaT", m_deltaT);
        }
    }

    // need to include verbose
    ParmParse ppAmr("AmrHydro");
    ppAmr.query("verbosity", m_verbosity);

    // BISICLES coupling
    ppParams.query("coupled_to_bisicles", m_coupled_to_bisicles);
    if (m_coupled_to_bisicles) {
        ppParams.get("bisicles_input_file", m_bisicles_input_file);
        ppParams.get("output_N_file", m_output_N_file);
        ppParams.query("ice_thickness_name", m_ice_thickness_name);
        ppParams.query("bed_elevation_name", m_bed_elevation_name);
        ppParams.query("velocity_x_name", m_velocity_x_name);
        ppParams.query("velocity_y_name", m_velocity_y_name);
    }

}

#include "NamespaceFooter.H"
