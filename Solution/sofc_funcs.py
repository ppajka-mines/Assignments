from math import exp
import numpy as np

"========= DEFINE RESIDUAL FUNCTION ========="
def residual(_, SV, pars, ptr):

    # Initialize the derivative / residual:
    dSV_dt = np.zeros_like(SV)

    # Anode potential does not change.  Leave those at zero.  Assumes very high
    #   electrical conductivity.


    # Anode double layer
    #   Overpotential:
    eta = -SV[ptr.phi_elyte[0]] - pars.E_an # Note phi_an = 0
    #   Faradaic current:
    i_Far_an = pars.i_o_an*(exp(-pars.aF_RT_an_fwd*eta) - exp(pars.aF_RT_an_rev*eta))
    #   Double-layer current:
    i_dl_an = -(pars.i_ext + i_Far_an)
    # This is for the electrolyte potential, which is equal to -dPhi_dl_an:
    dSV_dt[ptr.phi_elyte[0]] =  i_dl_an / pars.C_dl_an

    # Cathode double layer:
    #   Overpotential:
    eta = (SV[ptr.phi_ca[0]]-SV[ptr.phi_elyte[-1]]) - pars.E_ca
    #   Faradaic current:
    i_Far_ca = pars.i_o_ca*(exp(-pars.aF_RT_ca_fwd*eta) - exp(pars.aF_RT_ca_rev*eta))

    #   Double-layer current:
    i_dl_ca = pars.i_ext - i_Far_ca
    #   Cathode potential evolves at the rate of the local elyte, minus double layer:
    dSV_dt[ptr.phi_ca] =  - i_dl_ca / pars.C_dl_ca

    return dSV_dt