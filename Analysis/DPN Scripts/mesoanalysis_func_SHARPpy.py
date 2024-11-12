# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:38:07 2021

@author: wangm/dylan
"""

import numpy as np
import os
import sys
from metpy.units import units
# from metpy.plots import SkewT
import metpy.calc as mpcalc
from scipy.interpolate import interp1d

import sharppy
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo
"""
sfcpcl = params.parcelx( prof, flag=1 ) # Surface Parcel
fcstpcl = params.parcelx( prof, flag=2 ) # Forecast Parcel
mupcl = params.parcelx( prof, flag=3 ) # Most-Unstable Parcel
mlpcl = params.parcelx( prof, flag=4 ) # 100 mb Mean Layer Parcel
"""
def all_cape_cin(profset, sfcpcl, mlpcl, mupcl):
   
    #print('--------------SharpPy definition--------------')

    # Basic CAPE/CIN calculations using SHARPpy
    sbcape = sfcpcl.bplus * (units.joule/units.kilogram)
    sbcin = sfcpcl.bminus * (units.joule/units.kilogram)
    
    mlcape = mlpcl.bplus * (units.joule/units.kilogram)
    mlcin = mlpcl.bminus * (units.joule/units.kilogram)
    
    mucape = mupcl.bplus * (units.joule/units.kilogram)
    mucin = mupcl.bminus * (units.joule/units.kilogram)
    
    return sbcape, sbcin, mlcape, mlcin, mucape, mucin

def all_bkshears(profset):
    
    # For wind shears (Returned in kts according to doc), need to find the pressure at various level bounds
    sfcpres = profset.pres[profset.sfc]
    pres_1kmAGL = interp.pres(profset, interp.to_msl(profset, 1000.))
    pres_3kmAGL = interp.pres(profset, interp.to_msl(profset, 3000.))
    pres_6kmAGL = interp.pres(profset, interp.to_msl(profset, 6000.))
    
    # 0-1 km bulk shear. 
    bksh_sfc1km_U, bksh_sfc1km_V = winds.wind_shear(profset, pbot=sfcpres, ptop=pres_1kmAGL)
    bksh_sfc1km_dir = mpcalc.wind_direction(bksh_sfc1km_U*units.kt, bksh_sfc1km_V*units.kt)
    bksh_sfc1km_spd = mpcalc.wind_speed(bksh_sfc1km_U*units.kt, bksh_sfc1km_V*units.kt)
    
    # 0-3 km bulk shear
    bksh_sfc3km_U, bksh_sfc3km_V = winds.wind_shear(profset, pbot=sfcpres, ptop=pres_3kmAGL)
    bksh_sfc3km_dir = mpcalc.wind_direction(bksh_sfc3km_U*units.kt, bksh_sfc3km_V*units.kt)
    bksh_sfc3km_spd = mpcalc.wind_speed(bksh_sfc3km_U*units.kt, bksh_sfc3km_V*units.kt)
    
    # 0-6 km bulk shear
    bksh_sfc6km_U, bksh_sfc6km_V = winds.wind_shear(profset, pbot=sfcpres, ptop=pres_6kmAGL)
    bksh_sfc6km_dir = mpcalc.wind_direction(bksh_sfc6km_U*units.kt, bksh_sfc6km_V*units.kt)
    bksh_sfc6km_spd = mpcalc.wind_speed(bksh_sfc6km_U*units.kt, bksh_sfc6km_V*units.kt)
    
    # 1-6 km bulk shear
    bksh_1_6km_U, bksh_1_6km_V = winds.wind_shear(profset, pbot=pres_1kmAGL, ptop=pres_6kmAGL)
    bksh_1_6km_dir = mpcalc.wind_direction(bksh_1_6km_U*units.kt, bksh_1_6km_V*units.kt)
    bksh_1_6km_spd = mpcalc.wind_speed(bksh_1_6km_U*units.kt, bksh_1_6km_V*units.kt)
    
    # 1-3 km bulk shear
    bksh_1_3km_U, bksh_1_3km_V = winds.wind_shear(profset, pbot=pres_1kmAGL, ptop=pres_3kmAGL)
    bksh_1_3km_dir = mpcalc.wind_direction(bksh_1_3km_U*units.kt, bksh_1_3km_V*units.kt)
    bksh_1_3km_spd = mpcalc.wind_speed(bksh_1_3km_U*units.kt, bksh_1_3km_V*units.kt)
    
    # Effective-layer (for shear)
    # Get the most-unstable parcel EL height. Use a new most-unstable parcel definition instead
    # mupcl300mb_vals = params.DefineParcel(profset, flag=3, pres=300)     # the definition of pres depends on flag
    # mu_pcl300mb = params.cape(profset, lplvals=mu_pcl300mb_vals)       # This returns a parcel object, not cape
    # mupcl300mb = params.parcelx(profset, lplvals=mupcl300mb_vals)        # This works to give valid numbers for LCL, LFC, EL
    mupcl300mb = params.parcelx(profset, flag=3)
    
    parcel_topp = interp.pres(profset, interp.to_msl(profset, 500.))
    parcel_depth = sfcpres - parcel_topp
    
    mlpcl500m_vals = params.DefineParcel(profset, flag=4, pres=parcel_depth)
    mlpcl500m = params.parcelx(profset, lplvals=mlpcl500m_vals)
    
    # Calculate the effective inflow layer
    efflay_botpsh, efflay_toppsh = params.effective_inflow_layer(profset)
    
    #print(efflay_botpsh)
    if np.ma.isMaskedArray(efflay_botpsh) == False:
        # Effective layer exists, compute wind shear
        # Get the height AGL of the efflective layer base
        efflay_both = interp.to_agl(profset, interp.hght(profset, efflay_botpsh))
        efflay_depth = (mlpcl500m.elhght - efflay_both)/2
        # efflay_depth = (mupcl300mb.elhght - efflay_both)/2
                      # interp.pres(profset, interp.to_msl(profset, 6000.))
        efflay_toppsh = interp.pres(profset, interp.to_msl(profset, efflay_both+efflay_depth))
        bksh_eff_U, bksh_eff_V = winds.wind_shear(profset, pbot=efflay_botpsh, ptop=efflay_toppsh)           
        if np.ma.isMaskedArray(bksh_eff_U) == False:
            bksh_eff_dir = mpcalc.wind_direction(bksh_eff_U*units.kt, bksh_eff_V*units.kt)
            bksh_eff_spd = mpcalc.wind_speed(bksh_eff_U*units.kt, bksh_eff_V*units.kt)
        else:
            bksh_eff_dir = np.nan * units.deg
            bksh_eff_spd = np.nan * units.kt
    else:
        bksh_eff_U = np.nan
        bksh_eff_V = np.nan   
        bksh_eff_dir = np.nan * units.deg
        bksh_eff_spd = np.nan * units.kt
        
    return bksh_sfc1km_spd, bksh_sfc1km_dir, bksh_sfc3km_spd, bksh_sfc3km_dir, bksh_sfc6km_spd, bksh_sfc6km_dir, \
           bksh_eff_spd, bksh_eff_dir, bksh_1_6km_spd, bksh_1_6km_dir, bksh_1_3km_spd, bksh_1_3km_dir
           
def all_stormrel_helicity(profset):
    
    #print('------------SharpPy definition-----------')

    # Non-parcel-based storm motion experimental - output in kts
    #stMot_uR, stMot_vR, stMot_uL, stMot_vL = winds.non_parcel_bunkers_motion_experimental(profset)
    stMot_uR, stMot_vR, stMot_uL, stMot_vL = params.bunkers_storm_motion(profset)
    # print(stMot_uR, stMot_vR)

    # Mean wind calculations
    sfcpres = profset.pres[profset.sfc]
    pres_6kmAGL = interp.pres(profset, interp.to_msl(profset, 6000.))

    meanU, meanV = winds.mean_wind(profset, sfcpres, pres_6kmAGL)

    meanWindSpd = mpcalc.wind_speed(meanU*units.kt, meanV*units.kt)
    meanWindDir = mpcalc.wind_direction(meanU*units.kt, meanV*units.kt) 

    # For helicity or storm-relative helicity, use height above ground
    
    # Helicity - storm motion input in kts. Helicity output in m2/s2. bottom and top height
    # given in m AGL.
    SRH500m = winds.helicity(profset, 0, 500., stu=stMot_uR, stv=stMot_vR)[0] *  ((units.meter*units.meter)/(units.second*units.second))
    SRH1km = winds.helicity(profset, 0, 1000., stu=stMot_uR, stv=stMot_vR)[0] * ((units.meter*units.meter)/(units.second*units.second))
    SRH3km = winds.helicity(profset, 0, 3000., stu=stMot_uR, stv=stMot_vR)[0] * ((units.meter*units.meter)/(units.second*units.second))
    
    H500m = winds.helicity(profset, 0, 500., stu=0, stv=0)[0] * ((units.meter*units.meter)/(units.second*units.second))
    H1km = winds.helicity(profset, 0, 1000., stu=0, stv=0)[0] * ((units.meter*units.meter)/(units.second*units.second))
    H3km = winds.helicity(profset, 0, 3000., stu=0, stv=0)[0] * ((units.meter*units.meter)/(units.second*units.second))
    
    # Find the regular effective inflow layer (for SR helicity)
    efflay_botphel, efflay_topphel = params.effective_inflow_layer(profset)

    if np.ma.isMaskedArray(efflay_botphel) == False:
        efflay_bothhel = interp.to_agl(profset, interp.hght(profset, efflay_botphel))  # meters AGL
        efflay_tophhel = interp.to_agl(profset, interp.hght(profset, efflay_topphel))  # meters AGL
        ESRH = winds.helicity(profset, efflay_bothhel, efflay_tophhel, stu=stMot_uR, stv=stMot_vR)[0] * ((units.meter*units.meter)/(units.second*units.second))
        EH = winds.helicity(profset, efflay_bothhel, efflay_tophhel, stu=0, stv=0)[0] * ((units.meter*units.meter)/(units.second*units.second))
    else:
        efflay_bothhel = np.nan * units.meter
        efflay_tophhel = np.nan * units.meter
        ESRH = np.nan * ((units.meter*units.meter)/(units.second*units.second))
        EH = np.nan * ((units.meter*units.meter)/(units.second*units.second))
    
    stMotRM_Spd = mpcalc.wind_speed(stMot_uR*units.kt, stMot_vR*units.kt)      # in kts
    stMotRM_Dir = mpcalc.wind_direction(stMot_uR*units.kt, stMot_vR*units.kt)  # in deg
    
    efflay_basehgt = efflay_bothhel * units.meter
    efflay_tophgt = efflay_tophhel * units.meter
    efflay_depth = efflay_tophgt - efflay_basehgt
    
    return SRH500m, SRH1km, SRH3km, ESRH, H500m, H1km, H3km, EH, \
           stMotRM_Spd, stMotRM_Dir, meanWindSpd, meanWindDir, \
           efflay_basehgt, efflay_tophgt, efflay_depth

def all_lcllfcel(profset, sfcpcl, mlpcl, mupcl):
   
    #print('------------SharpPy definition---------------')

    # -- LCL, LFC, and EL. All values returned in m AGL.
    
    sb_LCL = sfcpcl.lclhght
    sb_LFC = sfcpcl.lfchght
    sb_EL = sfcpcl.elhght
    
    ml_LCL = mlpcl.lclhght
    ml_LFC = mlpcl.lfchght 
    ml_EL = mlpcl.elhght
    
    mu_LCL = mupcl.lclhght
    mu_LFC = mupcl.lfchght
    mu_EL = mupcl.elhght
    
    if np.ma.isMaskedArray(sb_LCL) == False:
        sb_LCL = sb_LCL * units.meter
    else:
        sb_LCL = np.nan * units.meter
    if np.ma.isMaskedArray(sb_LFC) == False:
        sb_LFC = sb_LFC * units.meter
    else:
        sb_LFC = np.nan * units.meter
    if np.ma.isMaskedArray(sb_EL) == False:
        sb_EL = sb_EL * units.meter
    else:
        sb_EL = np.nan * units.meter
    # -- 
    if np.ma.isMaskedArray(ml_LCL) == False:
        ml_LCL = ml_LCL * units.meter
    else:
        ml_LCL = np.nan * units.meter
    if np.ma.isMaskedArray(ml_LFC) == False:
        ml_LFC = ml_LFC * units.meter
    else:
        ml_LFC = np.nan * units.meter
    if np.ma.isMaskedArray(ml_EL) == False:
        ml_EL = ml_EL * units.meter
    else:
        ml_EL = np.nan * units.meter 
    # --
    if np.ma.isMaskedArray(mu_LCL) == False:
        mu_LCL = mu_LCL * units.meter
    else:
        mu_LCL = np.nan * units.meter
    if np.ma.isMaskedArray(mu_LFC) == False:
        mu_LFC = mu_LFC * units.meter
    else:
        mu_LFC = np.nan * units.meter
    if np.ma.isMaskedArray(mu_EL) == False:
        mu_EL = mu_EL * units.meter
    else:
        mu_EL = np.nan * units.meter    
        
        
    return sb_LCL, sb_LFC, sb_EL, ml_LCL, ml_LFC, ml_EL, mu_LCL, mu_LFC, mu_EL

def all_stormparams(profset, sfcpcl, mlpcl, mupcl):
    
    sbcape = all_cape_cin(profset, sfcpcl, mlpcl, mupcl)[0]
    sbcin = all_cape_cin(profset, sfcpcl, mlpcl, mupcl)[1]
    mlcape = all_cape_cin(profset, sfcpcl, mlpcl, mupcl)[2]
    mlcin = all_cape_cin(profset, sfcpcl, mlpcl, mupcl)[3]
    mucape = all_cape_cin(profset, sfcpcl, mlpcl, mupcl)[4]
    mucin = all_cape_cin(profset, sfcpcl, mlpcl, mupcl)[5]
    
    srh1km = all_stormrel_helicity(profset)[1]
    esrh = all_stormrel_helicity(profset)[3]
    
    sblcl = all_lcllfcel(profset, sfcpcl, mlpcl, mupcl)[0]
    mllcl = all_lcllfcel(profset, sfcpcl, mlpcl, mupcl)[3]
    
    bksh6km = all_bkshears(profset)[4]
    bksheff = all_bkshears(profset)[6]

    # check units of all input before calculating
    if sbcape.units != (units.joule/units.kilogram) or \
       sbcin.units != (units.joule/units.kilogram)  or \
       mlcape.units != (units.joule/units.kilogram) or \
       mlcin.units != (units.joule/units.kilogram) or \
       mucape.units != (units.joule/units.kilogram) or \
       mucin.units != (units.joule/units.kilogram):
       print('Warning: One of CAPE/CIN is not in the correct units!')
       sys.exit()
       
    if bksheff.units != units('m/s'):
        bksheff = bksheff.to(units('m/s'))
    if bksh6km.units != units('m/s'):
        bksh6km = bksh6km.to(units('m/s'))
    
    if esrh.units != ((units.meter*units.meter)/(units.second*units.second)) or \
       srh1km.units != ((units.meter*units.meter)/(units.second*units.second)):
       print('Warning: One of helicities is not in the correct units!')
       sys.exit()
       
    if sblcl.units != units.meter or \
       mllcl.units != units.meter:
       print('Warning: One of sounding levels is not in the correct units!')
       sys.exit()
    
    # Supercell composite parameter. SHARPpy is missing a CIN term. Adding it
    # myself.
    if mucin.m > -40:
        cin_term = 1
    else:
        cin_term = (-40/mucin.m)

    SCP_cin = params.scp(mucape.m, esrh.m, bksheff.m) * cin_term
    
    SCP = params.scp(mucape.m, esrh.m, bksheff.m)

    # Significant tornado parameter effective layer
    STP_eff = params.stp_cin(mlcape.m, esrh.m, bksheff.m, mllcl.m, mlcin.m)
    
    # Significant tornado parameter fixed-layer
    if sbcin.m > -50:
        cin_term = 1
    elif sbcin.m < -200:
        cin_term = 0
    else:
        cin_term = (200+sbcin.m)/150
    STP_fixed = params.stp_fixed(sbcape.m, sblcl.m, srh1km.m, bksh6km.m)*cin_term
    
    # Violent tornado parameter
    # I believe the 0-3 km CAPE is actually the MLCAPE integrated up to 3 km
    mlcape3km = mlpcl.b3km   # returned as J/kg
    lapse_rate3km = params.lapse_rate(profset, 0, 3000., pres=False)  # returned as degC/km, but not assigning a units here
    #if mlcape3km > 100:
    #    lapse_term = 2
    #else:
    #    lapse_term = lapse_rate3km/6.5
    
    lapse_term = lapse_rate3km/6.5

    mlcape3km_term = (mlcape3km/50.)
    if mlcape3km_term > 2:
        mlcape3km_term = 2
    VTP = STP_eff*mlcape3km_term*lapse_term
    
    #SHiP calculation for hail
    ship = params.ship(profset)
    
    return SCP, SCP_cin, STP_fixed, STP_eff, VTP, mlcape3km, lapse_rate3km, ship
    
def cape_10_40_calc(prof, pcl):
    cape = pcl.bplus * (units.joule/units.kilogram)
    cin = pcl.bminus * (units.joule/units.kilogram)
        
    return cape, cin