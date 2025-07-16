"""

Copyright 2019 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission
subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing,
software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.

"""
from __future__ import print_function, absolute_import
from nine import range

import warnings
import importlib
# importlib to import pandas as pd and excel for additional reservoirs as Excel table

from pcraster.operations import ifthen, boolean, defined, lookupscalar
import numpy as np


from ..global_modules.settings import LisSettings, MaskInfo
from ..global_modules.add1 import loadmap, compressArray, decompress, makenumpy, readwaterbody_Excel
from ..global_modules.errors import LisfloodWarning
from . import HydroModule




class reservoir(HydroModule):

    """
    # ************************************************************
    # ***** RESERVOIR    *****************************************
    # ************************************************************
    """
    input_files_keys = {'simulateReservoirs': ['ReservoirSites', 'TabTotStorage', 'TabConservativeStorageLimit',
                                              'TabNormalStorageLimit', 'TabFloodStorageLimit', 'TabNonDamagingOutflowQ',
                                              'TabNormalOutflowQ', 'TabMinOutflowQ', 'adjust_Normal_Flood',
                                              'ReservoirRnormqMult', 'ReservoirInitialFillValue']}
    module_name = 'Reservoir'

    def __init__(self, reservoir_variable):
        self.var = reservoir_variable

    def reservoir_releases(self,xl_settings_file_path):
        """
        Read daily reservoir release data as 2D numpy array [reservoirs ,366 days]
        :param self: all variables of self
        :param xl_settings_file_path: Excel file with sheet "Reservoirs_downstream"
        :return: 2D numpy array with release data per reservoir
        """

        pd = importlib.import_module ("pandas", package=None)
        df = pd.read_excel(xl_settings_file_path, sheet_name= 'Reservoirs_downstream')
        ReservoirSites_tolist = self.var.ReservoirSitesCC.tolist()

        reservoir_release = [[-1 for i in self.var.ReservoirSitesCC] for i in range(366)]
        for res in list(df)[2:]:
            if res in ReservoirSites_tolist:
                res_index = ReservoirSites_tolist.index(int(float(res)))

                for day in range(366):
                    reservoir_release[day][res_index] = df[res][day]

        return reservoir_release

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

    def initial(self):
        """ initial part of the reservoir module
        """
        # ************************************************************
        # ***** RESERVOIRS
        # ************************************************************
        settings = LisSettings.instance()
        option = settings.options
        maskinfo = MaskInfo.instance()
        if option['simulateReservoirs'] and not option['InitLisflood']:

            # NoSubStepsRes=max(1,roundup(self.var.DtSec/loadmap('DtSecReservoirs')))
            # Number of sub-steps based on value of DtSecReservoirs,
            # or 1 if DtSec is smaller than DtSecReservoirs
            # DtSubRes=self.var.DtSec/loadmap('NoSubStepsRes')
            # Corresponding sub-timestep (seconds)
            binding = settings.binding
            self.var.ReservoirSitesC = loadmap('ReservoirSites')
            self.var.ReservoirSitesC[self.var.ReservoirSitesC < 1] = 0
            self.var.ReservoirSitesC[self.var.IsChannel == 0] = 0

            # Get rid of any reservoirs that are not part of the channel network
            self.var.ReservoirSitesCC = np.compress(self.var.ReservoirSitesC > 0, self.var.ReservoirSitesC)
            self.var.ReservoirIndex = np.nonzero(self.var.ReservoirSitesC)[0]


            if option['reservoir_lakes_Excel']:
                xl_settings_file_path = "P:/watmodel/Lisflood/regions/morava_3arc/input/Lakes_Reservoirs_ID_tables/cwatm_settings_reservoirs_morava.xlsx"
                # look into the Excel to find if there are new reservoirs
                self.var.ReservoirSitesC = readwaterbody_Excel(self,xl_settings_file_path,self.var.ReservoirSitesC,2,3)

                self.var.ReservoirSitesCC = np.compress(self.var.ReservoirSitesC > 0, self.var.ReservoirSitesC)
                self.var.ReservoirIndex = np.nonzero(self.var.ReservoirSitesC)[0]


            # check if no reservoir is in the mask map area => no reservoir function is used
            if self.var.ReservoirSitesCC.size == 0:
                # break if no reservoirs
                warnings.warn(LisfloodWarning('There are no reservoirs. Reservoirs simulation won\'t run'))
                option['simulateReservoirs'] = False
                option['repsimulateReservoirs'] = False
                # rebuild lists of reported files with simulateReservoirs and repsimulateReservoirs = False
                settings.build_reportedmaps_dicts()
                return


            self.var.IsStructureKinematic = np.where(self.var.ReservoirSitesC > 0, np.bool8(1), self.var.IsStructureKinematic)
            # Add reservoir locations to structures map (used to modify LddKinematic
            # and to calculate LddStructuresKinematic)

            #ReservoirSitePcr = loadmap('ReservoirSites', pcr=True)
            ReservoirSitePcr = decompress(self.var.ReservoirSitesC)
            self.var.ReservoirSites = ReservoirSitePcr
            ReservoirSitePcr = ifthen((defined(ReservoirSitePcr) & boolean(decompress(self.var.IsChannel))), ReservoirSitePcr)
            # Get rid of any reservoirs that are not part of the channel network
            # (following logic of 'old' code the inflow into these reservoirs is
            # always zero, so either change this or leave them out!)

            TotalReservoirStorageM3 = lookupscalar(str(binding['TabTotStorage']), ReservoirSitePcr)
            self.var.TotalReservoirStorageM3C = compressArray(TotalReservoirStorageM3)
            self.var.TotalReservoirStorageM3C = np.where(np.isnan(self.var.TotalReservoirStorageM3C), 0, self.var.TotalReservoirStorageM3C)
            self.var.TotalReservoirStorageM3CC = np.compress(self.var.ReservoirSitesC > 0, self.var.TotalReservoirStorageM3C)
            # Total storage of each reservoir [m3]

            ConservativeStorageLimit = lookupscalar(str(binding['TabConservativeStorageLimit']), ReservoirSitePcr)
            ConservativeStorageLimitC = compressArray(ConservativeStorageLimit)
            self.var.ConservativeStorageLimitCC = np.compress(self.var.ReservoirSitesC > 0, ConservativeStorageLimitC)
            # Conservative storage limit (fraction of total storage, [-])

            NormalStorageLimit = lookupscalar(str(binding['TabNormalStorageLimit']), ReservoirSitePcr)
            NormalStorageLimitC = compressArray(NormalStorageLimit)
            self.var.NormalStorageLimitCC = np.compress(self.var.ReservoirSitesC > 0, NormalStorageLimitC)
            # Normal storage limit (fraction of total storage, [-])

            FloodStorageLimit = lookupscalar(str(binding['TabFloodStorageLimit']), ReservoirSitePcr)
            FloodStorageLimitC = compressArray(FloodStorageLimit)
            self.var.FloodStorageLimitCC = np.compress(self.var.ReservoirSitesC > 0, FloodStorageLimitC)
            # Flood storage limit (fraction of total storage, [-])

            NonDamagingReservoirOutflow = lookupscalar(str(binding['TabNonDamagingOutflowQ']), ReservoirSitePcr)
            NonDamagingReservoirOutflowC = compressArray(NonDamagingReservoirOutflow)
            self.var.NonDamagingReservoirOutflowCC = np.compress(self.var.ReservoirSitesC > 0, NonDamagingReservoirOutflowC)
            # Non-damaging reservoir outflow [m3/s]

            NormalReservoirOutflow = lookupscalar(str(binding['TabNormalOutflowQ']), ReservoirSitePcr)
            NormalReservoirOutflowC = compressArray(NormalReservoirOutflow)
            self.var.NormalReservoirOutflowCC = np.compress(self.var.ReservoirSitesC > 0, NormalReservoirOutflowC)
            # Normal reservoir outflow [m3/s]

            MinReservoirOutflow = lookupscalar(str(binding['TabMinOutflowQ']), ReservoirSitePcr)
            MinReservoirOutflowC = compressArray(MinReservoirOutflow)
            self.var.MinReservoirOutflowCC = np.compress(self.var.ReservoirSitesC > 0, MinReservoirOutflowC)
            # minimum reservoir outflow [m3/s]

            # Calibration
            adjust_Normal_Flood = loadmap('adjust_Normal_Flood')
            adjust_Normal_FloodC = makenumpy(adjust_Normal_Flood)
            adjust_Normal_FloodCC = np.compress(self.var.ReservoirSitesC > 0, adjust_Normal_FloodC)
            # adjusting the balance between normal and flood storage
            # big value (= closer to flood) results in keeping the normal qoutflow longer constant
            self.var.Normal_FloodStorageLimitCC = self.var.NormalStorageLimitCC + adjust_Normal_FloodCC * (self.var.FloodStorageLimitCC - self.var.NormalStorageLimitCC)

            ReservoirRnormqMult = loadmap('ReservoirRnormqMult')
            ReservoirRnormqMultC = makenumpy(ReservoirRnormqMult)
            ReservoirRnormqMultCC = np.compress(self.var.ReservoirSitesC > 0, ReservoirRnormqMultC)
            self.var.NormalReservoirOutflowCC = self.var.NormalReservoirOutflowCC * ReservoirRnormqMultCC
            # calibration: all reservoirs normal outflow are multiplied with a factor
            self.var.NormalReservoirOutflowCC = np.where(self.var.NormalReservoirOutflowCC > self.var.MinReservoirOutflowCC, self.var.NormalReservoirOutflowCC, self.var.MinReservoirOutflowCC+0.01)
            self.var.NormalReservoirOutflowCC = np.where(self.var.NormalReservoirOutflowCC < self.var.NonDamagingReservoirOutflowCC, self.var.NormalReservoirOutflowCC, self.var.NonDamagingReservoirOutflowCC-0.01)


            if option['reservoir_lakes_Excel']:
                # put here, because it comes before lakes and wetland initial
                self.var.waterbodyTypeCC = self.var.ReservoirSitesCC * 0

                for i in range(len(self.var.waterbody_info)):
                    resint = int(self.var.waterbody_info[i][0])
                    resindex = np.where(self.var.ReservoirSitesCC == resint)
                    # test if reservoir is found
                    if resindex[0].size > 0:
                        resindex = resindex[0].tolist()[0]
                        if int(self.var.waterbody_info[i][4]) > 0: self.var.waterbodyTypeCC[resindex] = int(self.var.waterbody_info[i][4])

                        if float(self.var.waterbody_info[i][7]) > 0: self.var.NormalReservoirOutflowCC[resindex] = float(self.var.waterbody_info[i][7])
                        # from Mio. m3 to m3 as rstor.txt is in m3
                        if float(self.var.waterbody_info[i][10]) > 0: self.var.TotalReservoirStorageM3CC[resindex] = float(self.var.waterbody_info[i][10]) * 1000000.

                        if float(self.var.waterbody_info[i][11]) > 0: self.var.ConservativeStorageLimitCC[resindex] = float(self.var.waterbody_info[i][11])
                        if float(self.var.waterbody_info[i][12]) > 0: self.var.NormalStorageLimitCC[resindex]  = float(self.var.waterbody_info[i][12])
                        if float(self.var.waterbody_info[i][13]) > 0: self.var.FloodStorageLimitCC[resindex]  = float(self.var.waterbody_info[i][13])
                        if float(self.var.waterbody_info[i][14]) > 0: self.var.MinReservoirOutflowCC[resindex] = float(self.var.waterbody_info[i][14])
                        if float(self.var.waterbody_info[i][15]) > 0: self.var.NormalReservoirOutflowCC[resindex] = float(self.var.waterbody_info[i][15])
                        if float(self.var.waterbody_info[i][16]) > 0: self.var.NonDamagingReservoirOutflowCC[resindex] = float(self.var.waterbody_info[i][16])

                        if float(self.var.waterbody_info[i][17]) > 0: adjust_Normal_FloodCC  = float(self.var.waterbody_info[i][17])
                        self.var.Normal_FloodStorageLimitCC[resindex]  = self.var.NormalStorageLimitCC[resindex]  + adjust_Normal_FloodCC * (
                                    self.var.FloodStorageLimitCC[resindex]  - self.var.NormalStorageLimitCC[resindex] )


            # Repeatedly used expressions in reservoirs routine
            self.var.DeltaO = self.var.NormalReservoirOutflowCC - self.var.MinReservoirOutflowCC
            self.var.DeltaLN = self.var.NormalStorageLimitCC - 2 * self.var.ConservativeStorageLimitCC
            self.var.DeltaLF = self.var.FloodStorageLimitCC - self.var.NormalStorageLimitCC
            self.var.DeltaNFL = self.var.FloodStorageLimitCC - self.var.Normal_FloodStorageLimitCC


            # read daily reservoir release data
            if option['reservoir_lakes_Excel'] & option['reservoir_release']:
                self.var.reservoir_releases = np.array(self.reservoir_releases(xl_settings_file_path))
                # self.var.CalendarDay

            
            ReservoirInitialFillValue = loadmap('ReservoirInitialFillValue')
            if np.max(ReservoirInitialFillValue) == -9999:
                ReservoirInitialFill = self.var.NormalStorageLimitCC,
            else:
                ReservoirInitialFill = np.compress(self.var.ReservoirSitesC > 0, ReservoirInitialFillValue)

            self.var.ReservoirFillCC = ReservoirInitialFill
            # Initial reservoir fill (fraction of total storage, [-])
            # -9999: assume reservoirs are filled to normal storage limit
            ReservoirStorageIniM3CC = ReservoirInitialFill * self.var.TotalReservoirStorageM3CC
            # Initial reservoir storage [m3] from state or initvalue
            self.var.ReservoirStorageM3CC = ReservoirStorageIniM3CC.copy()
            self.var.ReservoirFill = maskinfo.in_zero()
            # Initial fill of reservoirs (fraction of total storage, [-])

            self.var.ReservoirStorageIniM3 = maskinfo.in_zero()
            np.put(self.var.ReservoirStorageIniM3, self.var.ReservoirIndex, ReservoirStorageIniM3CC)

            self.var.ReservoirStorageM3 = self.var.ReservoirStorageIniM3


    def dynamic_inloop(self, NoRoutingExecuted):
        """ dynamic part of the lake routine
           inside the sub time step routing routine
        """

        # ************************************************************
        # ***** RESERVOIR
        # ************************************************************
        settings = LisSettings.instance()
        option = settings.options
        maskinfo = MaskInfo.instance()
        if option['simulateReservoirs'] and not option['InitLisflood']:

            InvDtSecDay = 1 / float(86400)
            # InvDtSecDay=self.var.InvDtSec
            # ReservoirInflow = cover(ifthen(defined(self.var.ReservoirSites), upstream(
            # self.var.LddStructuresKinematic, self.var.ChanQ)), scalar(0.0))

            ReservoirInflowCC = np.bincount(self.var.downstruct, weights=self.var.ChanQ)[self.var.ReservoirIndex]
            # ReservoirInflow=cover(ifpcr(defined(self.var.ReservoirSites),upstream(self.var.LddStructuresKinematic,self.var.ChanQ)),null)
            # Reservoir inflow in [m3/s]
            # 20-2-2006: Replaced ChanQKin by ChanQ (if this results in problems change back to ChanQKin!)
            # 21-2-2006: Inflow now taken from 1st upstream cell(s), using LddStructuresKinematic
            # (LddStructuresKinematic equals LddKinematic, but without the pits/sinks upstream of the structure
            # locations; note that using Ldd here instead would introduce MV!)

            QResInM3Dt = ReservoirInflowCC * self.var.DtRouting
            # Reservoir inflow in [m3] per timestep (routing step)
            
            # print('RESERVOIRS MODULE IN')
            # print(np.sum(self.var.ReservoirStorageM3))
            if NoRoutingExecuted==0:
                self.var.ReservoirStorageM3CC=np.compress(self.var.ReservoirSitesC > 0, self.var.ReservoirStorageM3)  ########################
            # print(np.sum(self.var.ReservoirStorageM3CC))
            
            self.var.ReservoirStorageM3CC += QResInM3Dt
            # New reservoir storage [m3] = plus inflow for this sub step
            self.var.ReservoirFillCC = self.var.ReservoirStorageM3CC / self.var.TotalReservoirStorageM3CC
            # New reservoir fill (fraction)

            ReservoirOutflow1 = np.minimum(self.var.MinReservoirOutflowCC, self.var.ReservoirStorageM3CC * InvDtSecDay)
            # Reservoir outflow [m3/s] if ReservoirFill le
            # 2*ConservativeStorageLimit

            ReservoirOutflow2 = self.var.MinReservoirOutflowCC + self.var.DeltaO * (self.var.ReservoirFillCC - 2 * self.var.ConservativeStorageLimitCC) / self.var.DeltaLN
            # Reservoir outflow [m3/s] if NormalStorageLimit le ReservoirFill
            # gt 2*ConservativeStorageLimit

            ReservoirOutflow3a = self.var.NormalReservoirOutflowCC
            ReservoirOutflow3b = self.var.NormalReservoirOutflowCC + ((self.var.ReservoirFillCC - self.var.Normal_FloodStorageLimitCC) / self.var.DeltaNFL) * (self.var.NonDamagingReservoirOutflowCC - self.var.NormalReservoirOutflowCC)
            # Reservoir outflow [m3/s] if FloodStorageLimit le ReservoirFill gt NormalStorageLimit
            # NEW 24-9-2004: linear transition between normal and non-damaging
            # outflow.
            #ReservoirOutflow4 = np.maximum((self.var.ReservoirFillCC - self.var.FloodStorageLimitCC) *
            #    self.var.TotalReservoirStorageM3CC * self.var.InvDtSec, self.var.NonDamagingReservoirOutflowCC)
            temp = np.minimum(self.var.NonDamagingReservoirOutflowCC, np.maximum(ReservoirInflowCC * 1.2, self.var.NormalReservoirOutflowCC))
            ReservoirOutflow4 = np.maximum((self.var.ReservoirFillCC - self.var.FloodStorageLimitCC-0.01) *
                                           self.var.TotalReservoirStorageM3CC * InvDtSecDay, temp)

            # Reservoir outflow [m3/s] if ReservoirFill gt FloodStorageLimit
            # Depending on ReservoirFill the reservoir outflow equals ReservoirOutflow1, ReservoirOutflow2,
            # ReservoirOutflow3 or ReservoirOutflow4

            ReservoirOutflow = ReservoirOutflow1.copy()
            ReservoirOutflow = np.where(self.var.ReservoirFillCC > 2 * self.var.ConservativeStorageLimitCC, ReservoirOutflow2, ReservoirOutflow)
           # ReservoirOutflow = np.where(self.var.ReservoirFillCC > self.var.NormalStorageLimitCC,
           #                     ReservoirOutflow3, ReservoirOutflow)

            ReservoirOutflow = np.where(self.var.ReservoirFillCC > self.var.NormalStorageLimitCC,
                                ReservoirOutflow3a, ReservoirOutflow)
            ReservoirOutflow = np.where(self.var.ReservoirFillCC > self.var.Normal_FloodStorageLimitCC,
                                ReservoirOutflow3b, ReservoirOutflow)

            ReservoirOutflow = np.where(self.var.ReservoirFillCC > self.var.FloodStorageLimitCC, ReservoirOutflow4, ReservoirOutflow)

            temp = np.minimum(ReservoirOutflow,np.maximum(ReservoirInflowCC, self.var.NormalReservoirOutflowCC))

            ReservoirOutflow = np.where((ReservoirOutflow > 1.2 * ReservoirInflowCC) &
                                        (ReservoirOutflow > self.var.NormalReservoirOutflowCC) &
                                        (self.var.ReservoirFillCC < self.var.FloodStorageLimitCC), temp, ReservoirOutflow)

            # new reservoir release function if reservoir watertyp is 3
            # if reservoir water typ is 2 (or not defined then use previous function
            if option['reservoir_lakes_Excel'] & option['reservoir_release']:
                reservoir_releaseCC = self.var.reservoir_releases[self.var.CalendarDay - 1]
                release = np.where(reservoir_releaseCC > -1, reservoir_releaseCC * self.var.ReservoirStorageM3CC * InvDtSecDay, ReservoirOutflow)
                ReservoirOutflow = np.where(self.var.ReservoirFillCC > self.var.FloodStorageLimitCC, ReservoirOutflow, release)

                #a = np.where(self.var.waterbody_typeCC == 3, 1.22, ReservoirOutflow)



            QResOutM3DtCC = ReservoirOutflow * self.var.DtRouting
            # Reservoir outflow in [m3] per sub step
            QResOutM3DtCC = np.minimum(QResOutM3DtCC, self.var.ReservoirStorageM3CC)
            # Check to prevent outflow from becoming larger than storage +
            # inflow
            QResOutM3DtCC = np.maximum(QResOutM3DtCC, self.var.ReservoirStorageM3CC - self.var.TotalReservoirStorageM3CC)

            # NEW 24-9-2004: Check to prevent reservoir storage from exceeding total capacity
            # expression to the right of comma always negative unless capacity is exceeded

            self.var.ReservoirStorageM3CC -= QResOutM3DtCC
            # New reservoir storage [m3]
            self.var.ReservoirFillCC = self.var.ReservoirStorageM3CC / self.var.TotalReservoirStorageM3CC
            # New reservoir fill

            # CM: Check ReservoirStorageM3CC for negative values and set them to zero
            self.var.ReservoirFillCC[np.isnan(self.var.ReservoirFillCC)] = 0
            self.var.ReservoirFillCC[self.var.ReservoirFillCC < 0] = 0

            # nel = len(self.var.ReservoirFillCC[:])  # always 1 # FIXME
            # for i in range(0, nel-1):  # never entering in loop
            #     if np.isnan(self.var.ReservoirFillCC[i]) or self.var.ReservoirFillCC[i] < 0:
            #         msg = "Negative or NaN volume for reservoir fill set to 0. Increase computation time step for routing (DtSecChannel) \n"
            #         warnings.warn(LisfloodWarning(msg))
            #         self.var.ReservoirFillCC[self.var.ReservoirFillCC < 0] = 0
            #         self.var.ReservoirFillCC[np.isnan(self.var.ReservoirFillCC)] = 0

            # Check ReservoirFillCC for negative values and set them to zero
            # if np.isnan(self.var.ReservoirFillCC).any() or (self.var.ReservoirFillCC < 0).any():
            #     msg = "Negative or NaN volume for reservoir fill set to 0. Increase computation time step for routing (DtSecChannel)"
            #     warnings.warn(LisfloodWarning(msg))
            #     self.var.ReservoirFillCC[self.var.ReservoirFillCC < 0] = 0
            #     self.var.ReservoirFillCC[np.isnan(self.var.ReservoirFillCC)] = 0

            # expanding the size as input for routing routine
            self.var.QResOutM3Dt = maskinfo.in_zero()
            np.put(self.var.QResOutM3Dt,self.var.ReservoirIndex,QResOutM3DtCC)
            # this is put to the channel again at each sub timestep

            if option['repsimulateReservoirs']:
                if NoRoutingExecuted == 0:
                    self.var.ReservoirInflowM3S = maskinfo.in_zero()
                    self.var.ReservoirOutflowM3S = maskinfo.in_zero()
                    self.var.sumResInCC = QResInM3Dt
                    self.var.sumResOutCC = QResOutM3DtCC
                    # for timeseries output - in and outflow to the reservoir is sumed up over the sub timesteps and stored in m/s
                    # set to zero at first timestep
                else:
                    self.var.sumResInCC += QResInM3Dt
                    self.var.sumResOutCC += QResOutM3DtCC
                    # summing up over all sub timesteps



            if NoRoutingExecuted == (self.var.NoRoutSteps-1):

                # expanding the size after last sub timestep
                self.var.ReservoirStorageM3 = maskinfo.in_zero()
                self.var.ReservoirFill = maskinfo.in_zero()
                np.put(self.var.ReservoirStorageM3, self.var.ReservoirIndex, self.var.ReservoirStorageM3CC)
                np.put(self.var.ReservoirFill, self.var.ReservoirIndex, self.var.ReservoirFillCC)

                # print('RESERVOIRS MODULE OUT')
                # print(np.sum(self.var.ReservoirStorageM3))
                # print(np.sum(self.var.ReservoirStorageM3CC))


                if option['repsimulateReservoirs']:
                    np.put(self.var.ReservoirInflowM3S, self.var.ReservoirIndex, self.var.sumResInCC / self.var.DtSec)
                    np.put(self.var.ReservoirOutflowM3S, self.var.ReservoirIndex, self.var.sumResOutCC / self.var.DtSec)
