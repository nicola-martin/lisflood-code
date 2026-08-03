"""

Copyright 2019 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""
from __future__ import absolute_import, print_function

import warnings
import importlib
# importlib to import pandas as pd and excel for additional reservoirs as Excel table
import datetime


import numpy as np
import pcraster

from ..global_modules.errors import LisfloodWarning, LisfloodError
from ..global_modules.add1 import loadmap, compressArray, decompress, readwaterbody_Excel
from ..global_modules.settings import LisSettings, MaskInfo
from . import HydroModule


class wetlands(HydroModule):

    """
    # ************************************************************
    # ***** WETLANDS     *****************************************
    # ************************************************************
    """

    input_files_keys = {'simulateWetlands': ['WetlandInitialLevelValue', 'TabWetlandAvNetInflowEstimate', 'PrevDischarge',
                                          'WetlandPrevInflowValue', 'WetlandPrevOutflowValue']}

    module_name = 'Wetlands'

    def __init__(self, wetlands_variable):
        self.var = wetlands_variable

# --------------------------------------------------------------------------
    def wetland_readarea(self, xl_settings_file_path):
        """
        Read daily wetland area as 2D numpy array [wetland area per wetland ,366 days]
        :param self: all variables of self
        :param xl_settings_file_path: Excel file with sheet "Wetlands"
        :return: 2D numpy array wetland area
        """

        pd = importlib.import_module("pandas", package=None)
        try:
            df = pd.read_excel(xl_settings_file_path, sheet_name='Wetlands')
        except:
            msg ="The Excel file: {} exists but it does not have the sheet: Wetlands\n".format(xl_settings_file_path)
            raise LisfloodError(msg)

        wetlandSites_tolist = self.var.WetlandSitesCC.tolist()

        # initialize wetlands for all 366 days
        wetland_area = np.zeros((366,len(self.var.WetlandSitesCC))) -1

        # Excel sheet from column 5 ->
        for res in list(df)[5:]:
            if res in wetlandSites_tolist:
                wet_index = wetlandSites_tolist.index(int(float(res)))
                for day in range(366):
                    wetland_area[day][wet_index] = df[res][day+3]

        return np.array(wetland_area)


# --------------------------------------------------------------------------

    def initial(self):
        """ initial part of the wetlands module
        """

        # ************************************************************
        # ***** WETLAND
        # ************************************************************
        settings = LisSettings.instance()
        option = settings.options
        binding = settings.binding
        maskinfo = MaskInfo.instance()

        # optional Wetland in settings. If not in settings do not use wetlands
        # First time it is called it sets a False if not in option
        if not('simulateWetlands' in option):
            option['simulateWetlands'] = False
        if option['simulateWetlands'] and not option['InitLisflood']:

            self.var.WetlandSitesC = loadmap('WetlandSites') * 0
            self.var.WetlandSitesC[self.var.WetlandSitesC < 1] = 0
            self.var.WetlandSitesC[self.var.IsChannel == 0] = 0
            # Get rid of any wetlands that are not part of the channel network

            # mask wetlands sites when using sub-catchments mask
            self.var.WetlandSitesCC = np.compress(self.var.WetlandSitesC > 0, self.var.WetlandSitesC)
            self.var.WetlandIndex = np.nonzero(self.var.WetlandSitesC)[0]

            # make it optional: if there is no option in the seetingsfile put it to False
            if not('reservoir_lakes_Excel' in option):
                option['reservoir_lakes_Excel'] = False
            if option['reservoir_lakes_Excel']:
                # load location of Excel file from settingsflie
                self.var.xl_settings_file_path = binding['Excel_settings_file']
                # look into the Excel to find if there are new wetland => waterbodyType = 5
                self.var.WetlandSitesC  = readwaterbody_Excel(self,self.var.xl_settings_file_path,self.var.WetlandSitesC, 5,5)
                self.var.WetlandSitesCC = np.compress(self.var.WetlandSitesC > 0, self.var.WetlandSitesC)
                self.var.WetlandIndex = np.nonzero(self.var.WetlandSitesC)[0]
            else:
                msg = "Wetland routine needs option: reservoir_lakes_Excel = True and sheet: Wetlands"
                raise LisfloodError(msg)


            if self.var.WetlandSitesCC.size == 0:
                #warnings.warn(LisfloodWarning('There are no wetlands. Wetlands simulation won\'t run'))
                option['simulateWetlands'] = False
                option['repsimulateWetlands'] = False
                # rebuild lists of reported files with repsimulateWetlands and simulateWetlands = False
                settings.build_reportedmaps_dicts()
                return
            # break if no wetlands

            self.var.IsStructureKinematic = np.where(self.var.WetlandSitesC > 0, np.bool8(1), self.var.IsStructureKinematic)
            # Add wetland locations to structures map (used to modify LddKinematic
            # and to calculate LddStructuresKinematic)

            # PCRaster part
            # -----------------------
            #WetlandSitePcr = loadmap('WetlandSites', pcr=True)
            # not need because it is decompressing the numpy arrea (which can have additional wetland from Excel file)
            WetlandSitePcr = decompress(self.var.WetlandSitesC)
            WetlandSitePcr = pcraster.ifthen((pcraster.defined(WetlandSitePcr) & pcraster.boolean(decompress(self.var.IsChannel))), WetlandSitePcr)
            IsStructureWetland = pcraster.boolean(WetlandSitePcr)
            # additional structure map only for wetlands to calculate water balance
            self.var.IsUpsOfStructureWetland = pcraster.downstream(self.var.LddKinematic, pcraster.cover(IsStructureWetland, 0))
            # Get all pixels just upstream of wetlands
            # -----------------------

            self.var.WetlandInflowOldCC = np.bincount(self.var.downstruct, weights=self.var.ChanQ)[self.var.WetlandIndex]
            # for Modified Puls Method the Q(inflow)1 has to be used.
            # It is assumed that this is the same as Q(inflow)2 for the first timestep
            # has to be checked if this works in forecasting mode!

            # -> can be changed, but now add the area to the lake table (if you use the ezxcel only this is not necessary)
            WetlandArea = pcraster.lookupscalar(str(binding['TabLakeArea']), WetlandSitePcr)
            WetlandAreaC = compressArray(WetlandArea)
            self.var.WetlandAreaCC = np.compress(self.var.WetlandSitesC > 0, WetlandAreaC)
                                 
            # Surface area of each wetland [m2]
            WetlandA = pcraster.lookupscalar(str(binding['TabLakeA']), WetlandSitePcr)
            # -> it is the same  for lakes(not an additional for wetlands only)
            WetlandMult = loadmap('LakeMultiplier')
            WetlandAC = compressArray(WetlandA) * WetlandMult
            self.var.WetlandACC = np.compress(self.var.WetlandSitesC > 0, WetlandAC)
            # Wetland parameter A (suggested  value equal to outflow width in [m])
            # multiplied with the calibration parameter WetlandMultiplier

            maskone = maskinfo.in_zero() + 1
            self.var.WetlandETmult  = np.compress(self.var.WetlandSitesC > 0, maskone * loadmap('WetlandETmult'))
            # wetland actual evapotranspiration multiplier
            self.var.wetland_maxlevel = np.compress(self.var.WetlandSitesC > 0, maskone * loadmap('Wetland_maxlevel'))


            if option['reservoir_lakes_Excel']:
                # if wetlands are stored in Excel file
                self.var.wetland_area = self.wetland_readarea(self.var.xl_settings_file_path)
                doy = self.var.CalendarDay.timetuple().tm_yday
                self.var.WetlandAreaCC = self.var.wetland_area[doy-1,:] * 1000000
                # back to wetlandArea , because it is used in routing_kinematic
                self.var.wetlandArea = maskinfo.in_zero()
                np.put(self.var.wetlandArea, self.var.WetlandIndex, self.var.WetlandAreaCC)

                for i in range(len(self.var.waterbody_info)):
                    wetlandint = int(self.var.waterbody_info[i][0])
                    wetlandindex = np.where(self.var.WetlandSitesCC == wetlandint)
                    # test if reservoir is found
                    if wetlandindex[0].size > 0:
                        wetlandindex = wetlandindex[0].tolist()[0]
                        # from km2 to m2 as wetlandarea.txt is in m2
                        ## outcommented because flexible area is used for wetlands
                        ###if float(self.var.waterbody_info[i][6]) > 0: self.var.WetlandAreaCC[wetlandindex] = float(self.var.waterbody_info[i][6]) * 1000000.
                        # if lake multiplier is in Excel use this -> it is the same for lakes (not an additional for wetlands only)
                        if float(self.var.waterbody_info[i][9]) > 0:
                            mult = float(self.var.waterbody_info[i][9])
                        else:
                            mult = np.compress(self.var.WetlandSitesCC > 0, maskone * WetlandMult)
                            mult = mult[wetlandindex]

                        # if wetland maxlevel is in Excel use this
                        if float(self.var.waterbody_info[i][18]) > 0:
                            self.var.wetland_maxlevel[wetlandindex]  = float(self.var.waterbody_info[i][18])
                        if float(self.var.waterbody_info[i][8]) > 0:
                            self.var.WetlandACC[wetlandindex]  = float(self.var.waterbody_info[i][8]) * mult
                        else:
                            # if no wetlandA is given then use Normal discharge to calculate LakeA
                            chanwidth = 7.1 * np.power(float(self.var.waterbody_info[i][7]), 0.539)
                            self.var.WetlandACC[wetlandindex] = mult * 0.612 * 2 / 3 * chanwidth * (2 * 9.81) ** 0.5


            WetlandInitialLevelValue  = loadmap('WetlandInitialLevelValue')
            if np.max(WetlandInitialLevelValue) == -9999:
                WetlandAvNetInflowEstimate = pcraster.lookupscalar(str(binding['TabWetlandAvNetInflowEstimate']), WetlandSitePcr)
                WetlandAvNetC = compressArray(WetlandAvNetInflowEstimate)
                self.var.WetlandAvNetCC = np.compress(self.var.WetlandSitesC > 0, WetlandAvNetC)

                WetlandStorageIniM3CC = self.var.WetlandAreaCC * np.sqrt(self.var.WetlandAvNetCC / self.var.WetlandACC)
                # Initial wetland storage [m3]  based on: S = WetlandArea * H = WetlandArea * sqrt(Q/a)
                self.var.WetlandLevelCC = WetlandStorageIniM3CC / self.var.WetlandAreaCC
            else:
                self.var.WetlandLevelCC = np.compress(self.var.WetlandSitesC > 0, WetlandInitialLevelValue)
                WetlandStorageIniM3CC = self.var.WetlandAreaCC * self.var.WetlandLevelCC
                # Initial wetland storage [m3]  based on: S = WetlandArea * H
                self.var.WetlandAvNetCC = np.compress(self.var.WetlandSitesC > 0, loadmap('PrevDischarge'))

            WetlandPrevInflowValue  = loadmap('WetlandPrevInflowValue')
            if np.max(WetlandInitialLevelValue) == -9999:
                self.var.WetlandInflowOldCC = np.bincount(self.var.downstruct, weights = self.var.ChanQ)[self.var.WetlandIndex]
            else:
                self.var.WetlandInflowOldCC = np.compress(self.var.WetlandSitesC > 0, WetlandPrevInflowValue)

            # Repeatedly used expressions in wetland routine

            # NEW Wetland Routine using Modified Puls Method (see Maniak, p.331ff)
            # (Qin1 + Qin2)/2 - (Qout1 + Qout2)/2 = (S2 - S1)/dtime
            # changed into:
            # (S2/dtime + Qout2/2) = (S1/dtime + Qout1/2) - Qout1 + (Qin1 + Qin2)/2
            # outgoing discharge (Qout) are linked to storage (S) by elevation.
            # Now some assumption to make life easier:
            # 1.) storage volume is increase proportional to elevation: S = A * H
            #      H: elevation, A: area of wetland
            # 2.) outgoing discharge = c * b * H **2.0 (c: weir constant, b: width)
            #      2.0 because it fits to a parabolic cross section see Aigner 2008
            #      (and it is much easier to calculate (that's the main reason)
            # c for a perfect weir with mu=0.577 and Poleni: 2/3 mu * sqrt(2*g) = 1.7
            # c for a parabolic weir: around 1.8
            # because it is a imperfect weir: C = c* 0.85 = 1.5
            # results in a formular : Q = 1.5 * b * H ** 2 = a*H**2 -> H =
            # sqrt(Q/a)
            self.var.WetlandFactor = self.var.WetlandAreaCC / (self.var.DtRouting * np.sqrt(self.var.WetlandACC))

            #  solving the equation  (S2/dtime + Qout2/2) = (S1/dtime + Qout1/2) - Qout1 + (Qin1 + Qin2)/2
            #  SI = (S2/dtime + Qout2/2) =  (A*H)/DtRouting + Q/2 = A/(DtRouting*sqrt(a)  * sqrt(Q) + Q/2
            #  -> replacement: A/(DtRouting*sqrt(a)) = Wetlandfactor, Y = sqrt(Q)
            #  Y**2 + 2*Wetlandfactor*Y-2*SI=0
            # solution of this quadratic equation:
            # Q=sqr(-WetlandFactor+sqrt(sqr(WetlandFactor)+2*SI))

            self.var.WetlandFactorSqr = np.square(self.var.WetlandFactor)
            # for faster calculation inside dynamic section
            WetlandStorageIndicator = WetlandStorageIniM3CC / self.var.DtRouting + self.var.WetlandAvNetCC / 2
            # SI = S/dt + Q/2

            WetlandPrevOutflowValue  = loadmap('WetlandPrevOutflowValue')
            if np.max(WetlandPrevOutflowValue) == -9999:
                # wetland lakes are not rectagular but have a triangular shape
                # Instead of lake -> wetland formular is a bit different:
                ##self.var.WetlandOutflowCC = np.square(-self.var.WetlandFactor + np.sqrt(self.var.WetlandFactorSqr + 2 * WetlandStorageIndicator))
                self.var.WetlandOutflowCC = np.square(-0.5 * self.var.WetlandFactor + np.sqrt(0.25 * self.var.WetlandFactorSqr + 2 * WetlandStorageIndicator))

                # solution of quadratic equation
                # it is as easy as this because:
                # 1. storage volume is increase proportional to elevation
                # 2. Q= a *H **2.0  (if you choose Q= a *H **1.5 you have to solve
                # the formula of Cardano)
            else:
                self.var.WetlandOutflowCC = np.compress(self.var.WetlandSitesC > 0, WetlandPrevOutflowValue)


            # if not initialized before, initialize excel wetlands
            if option['reservoir_lakes_Excel']:
                for i in range(len(self.var.waterbody_info)):
                    wetlandint = int(self.var.waterbody_info[i][0])
                    wetlandindex = np.where(self.var.WetlandSitesCC == wetlandint)
                    # test if reservoir is found
                    if wetlandindex[0].size > 0:
                        wetlandindex = wetlandindex[0].tolist()[0]
                        if (np.max(WetlandInitialLevelValue) == -9999) | (np.isnan(self.var.WetlandAvNetCC[wetlandindex])):
                            self.var.WetlandAvNetCC[wetlandindex] = float(self.var.waterbody_info[i][7])
                            WetlandStorageIniM3CC[wetlandindex] = self.var.WetlandAreaCC[wetlandindex] * np.sqrt(self.var.WetlandAvNetCC[wetlandindex] / self.var.WetlandACC[wetlandindex])
                            self.var.WetlandLevelCC[wetlandindex] = WetlandStorageIniM3CC[wetlandindex] / self.var.WetlandAreaCC[wetlandindex]
                            WetlandStorageIndicator[wetlandindex]  = WetlandStorageIniM3CC[wetlandindex]  / self.var.DtRouting + self.var.WetlandAvNetCC[wetlandindex]  / 2
                            self.var.WetlandOutflowCC[wetlandindex]  = np.square(-0.5* self.var.WetlandFactor[wetlandindex]  + np.sqrt(0.25 * self.var.WetlandFactorSqr[wetlandindex]  + 2 * WetlandStorageIndicator[wetlandindex]))

            self.var.WetlandStorageM3CC = WetlandStorageIniM3CC.copy()
            self.var.WetlandStorageM3BalanceCC = WetlandStorageIniM3CC.copy()

            # lake level is average lake level = 1/2 of  max level for a triangular lake
            #  lakelevel should be at wetland_maxlevel (e.g. =1.0 m) -> rest goes to outflow
            # if lakelevel >= maxlevel (default =1.0) sea level is kept constant and equation is changing
            wetlandOut2 = np.maximum(0, 2 * (WetlandStorageIndicator- self.var.wetland_maxlevel * self.var.WetlandAreaCC/self.var.DtRouting))
            self.var.WetlandOutflowCC = np.where((self.var.WetlandLevelCC >= self.var.wetland_maxlevel), wetlandOut2,self.var.WetlandOutflowCC)

            self.var.WetlandStorageIniM3 = maskinfo.in_zero()
            self.var.WetlandLevel = maskinfo.in_zero()
            self.var.WetlandInflowOld = maskinfo.in_zero()
            self.var.WetlandOutflow = maskinfo.in_zero()
            np.put(self.var.WetlandStorageIniM3,self.var.WetlandIndex,WetlandStorageIniM3CC)
            self.var.WetlandStorageM3 = self.var.WetlandStorageIniM3.copy()
            np.put(self.var.WetlandLevel, self.var.WetlandIndex, self.var.WetlandLevelCC)
            np.put(self.var.WetlandInflowOld, self.var.WetlandIndex, self.var.WetlandInflowOldCC)
            np.put(self.var.WetlandOutflow, self.var.WetlandIndex, self.var.WetlandOutflowCC)

            self.var.EWWetlandCUMM3 = maskinfo.in_zero()
            self.var.EWWetlandWBM3 = maskinfo.in_zero()
            # Initialising cumulative output variables
            # These are all needed to compute the cumulative mass balance error



    def dynamic_inloop(self, NoRoutingExecuted):
        """ dynamic part of the wetland routine
           inside the sub time step routing routine
        """

        # ************************************************************
        # ***** WETLAND
        # ************************************************************
        settings = LisSettings.instance()
        option = settings.options
        maskinfo = MaskInfo.instance()

        if option['simulateWetlands'] and not(option['InitLisflood']):    # only with no InitLisflood

            if NoRoutingExecuted==0:
                self.var.WetlandAreaCC = self.var.wetland_area[self.var.CalendarDay - 1, :] * 1000000
                np.put(self.var.wetlandArea, self.var.WetlandIndex, self.var.WetlandAreaCC)

                self.var.WetlandStorageM3CC=np.compress(self.var.WetlandSitesC > 0, self.var.WetlandStorageM3)
                # Evaporation from lakes is calculated
                ewWetlandCC = np.compress(self.var.WetlandSitesC > 0, self.var.EWRef)
                # higher evaporation due to waterplants and swallow lakes

                ewWetlandCC = ewWetlandCC * self.var.WetlandETmult

                # evaporation from open water [mm] to [m] * lake area [m2]
                self.var.evaWetlandCC = (ewWetlandCC * 0.001 * self.var.WetlandAreaCC) / self.var.NoRoutSteps
                self.var.evaWetlandCCsum = self.var.evaWetlandCC * 0

            self.var.WetlandInflowCC = np.bincount(self.var.downstruct, weights=self.var.ChanQ)[self.var.WetlandIndex]
            # Wetland inflow in [m3/s]

            WetlandIn = (self.var.WetlandInflowCC + self.var.WetlandInflowOldCC) * 0.5
            # for Modified Puls Method: (S2/dtime + Qout2/2) = (S1/dtime + Qout1/2) - Qout1 + (Qin1 + Qin2)/2
            #  here: (Qin1 + Qin2)/2
            self.var.WetlandInflowOldCC = self.var.WetlandInflowCC.copy()
            # Qin2 becomes Qin1 for the next time step

            # calculate real evaporation (not more then water in the lake)
            evaWetlandCC = np.where((self.var.WetlandStorageM3CC - self.var.evaWetlandCC) > 0., self.var.evaWetlandCC, self.var.WetlandStorageM3CC)
            # sum up real evaporation from lake
            self.var.evaWetlandCCsum += evaWetlandCC
            # Lake storage minus lake evaporation
            self.var.WetlandStorageM3CC = self.var.WetlandStorageM3CC - evaWetlandCC

            WetlandStorageIndicator = self.var.WetlandStorageM3CC /self.var.DtRouting - 0.5 * self.var.WetlandOutflowCC + WetlandIn
            # here S1/dtime - Qout1/2 + WetlandIn , so that is the right part
            # of the equation above

            # calculation if var.waterBodyTyp = 5 and lake is assumed to be triangular
            # and therefore the equation is a bit different
            #self.var.WetlandOutflowCC = np.square( -self.var.WetlandFactor + np.sqrt(self.var.WetlandFactorSqr + 2 * WetlandStorageIndicator))
            self.var.WetlandOutflowCC = np.square(-0.5 * self.var.WetlandFactor + np.sqrt(0.25 * self.var.WetlandFactorSqr + 2 * WetlandStorageIndicator))

            #  lakelevel should be at wetland_maxlevel (default =1.0 m) -> rest goes to outflow
            # if lakelevel >= maxlevel sea level is kept constant and equation is changing
            testlevel = ((WetlandStorageIndicator - self.var.WetlandOutflowCC * 0.5) * self.var.DtRouting) / self.var.WetlandAreaCC
            #outflow adjusted to reach self.var.wetland_maxlevel
            wetlandOut2 = np.maximum(0, 2 * (WetlandStorageIndicator- self.var.wetland_maxlevel * self.var.WetlandAreaCC/self.var.DtRouting))
            self.var.WetlandOutflowCC = np.where((testlevel > self.var.wetland_maxlevel), wetlandOut2,self.var.WetlandOutflowCC)

            # Flow out of wetland:
            #  solving the equation  (S2/dtime + Qout2/2) = (S1/dtime + Qout1/2) - Qout1 + (Qin1 + Qin2)/2
            #  SI = (S2/dtime + Qout2/2) =  (A*H)/DtRouting + Q/2 = A/(DtRouting*sqrt(a)  * sqrt(Q) + Q/2
            #  -> replacement: A/(DtRouting*sqrt(a)) = Wetlandfactor, Y = sqrt(Q)
            #  Y**2 + 2*Wetlandfactor*Y-2*SI=0
            # solution of this quadratic equation:
            # Q=sqr(-WetlandFactor+sqrt(sqr(WetlandFactor)+2*SI));

            # expanding the size to save as state variable
            self.var.WetlandOutflow = maskinfo.in_zero()
            np.put(self.var.WetlandOutflow, self.var.WetlandIndex, self.var.WetlandOutflowCC)

            QWetlandOutM3DtCC = self.var.WetlandOutflowCC * self.var.DtRouting
            # Outflow in [m3] per timestep
            # Needed at every cell, hence cover statement

            self.var.WetlandStorageM3CC = (WetlandStorageIndicator - self.var.WetlandOutflowCC* 0.5) * self.var.DtRouting
            # Wetland storage

            # self.var.WetlandStorageM3CC < 0 leads to NaN in state files
            # Check WetlandStorageM3CC for negative values and set them to zero
            if any(np.isnan(self.var.WetlandStorageM3CC)) or any(self.var.WetlandStorageM3CC < 0):
                msg = "Negative or NaN volume for wetland storage set to 0. " \
                      "Increase computation time step for routing (DtSecChannel) \n"
                warnings.warn(LisfloodWarning(msg))
                self.var.WetlandStorageM3CC[self.var.WetlandStorageM3CC < 0] = 0
                self.var.WetlandStorageM3CC[np.isnan(self.var.WetlandStorageM3CC)] = 0

            self.var.WetlandStorageM3BalanceCC += WetlandIn * self.var.DtRouting - QWetlandOutM3DtCC - evaWetlandCC
            # for mass balance, the wetland storage is calculated every time step
            self.var.WetlandLevelCC = self.var.WetlandStorageM3CC / self.var.WetlandAreaCC

            # expanding the size
            self.var.QWetlandOutM3Dt = maskinfo.in_zero()
            np.put(self.var.QWetlandOutM3Dt,self.var.WetlandIndex,QWetlandOutM3DtCC)

            if option['repsimulateWetlands']:
                if NoRoutingExecuted == 0:
                    self.var.WetlandInflowM3S = maskinfo.in_zero()
                    self.var.WetlandOutflowM3S = maskinfo.in_zero()
                    self.var.sumWetlandInCC = self.var.WetlandInflowCC * self.var.DtRouting
                    self.var.sumWetlandOutCC = QWetlandOutM3DtCC
                    # for timeseries output - in and outflow to the reservoir
                    # is sumed up over the sub timesteps and stored in m/s
                    # set to zero at first timestep
                else:
                    self.var.sumWetlandInCC += self.var.WetlandInflowCC * self.var.DtRouting
                    self.var.sumWetlandOutCC += QWetlandOutM3DtCC
                    # summing up over all sub timesteps

            if NoRoutingExecuted == (self.var.NoRoutSteps-1):
                # expanding the size after last sub timestep
                self.var.WetlandStorageM3Balance = maskinfo.in_zero()
                self.var.WetlandStorageM3 = maskinfo.in_zero()
                self.var.WetlandLevel = maskinfo.in_zero()
                self.var.WetlandInflowOld = maskinfo.in_zero()
                self.var.WetlandOutflow = maskinfo.in_zero()
                self.var.evaWetlandM3 = maskinfo.in_zero()
                np.put(self.var.WetlandStorageM3Balance, self.var.WetlandIndex, self.var.WetlandStorageM3BalanceCC)
                np.put(self.var.WetlandStorageM3, self.var.WetlandIndex, self.var.WetlandStorageM3CC)
                np.put(self.var.WetlandLevel, self.var.WetlandIndex, self.var.WetlandLevelCC)
                np.put(self.var.WetlandInflowOld, self.var.WetlandIndex, self.var.WetlandInflowOldCC)
                np.put(self.var.WetlandOutflow, self.var.WetlandIndex, self.var.WetlandOutflowCC)
                np.put(self.var.evaWetlandM3, self.var.WetlandIndex, self.var.evaWetlandCCsum)

                if option['repsimulateWetlands']:
                    np.put(self.var.WetlandInflowM3S, self.var.WetlandIndex, self.var.sumWetlandInCC / self.var.DtSec)
                    np.put(self.var.WetlandOutflowM3S, self.var.WetlandIndex, self.var.sumWetlandOutCC / self.var.DtSec)
