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
from nine import range

from pcraster import ifthenelse, downstream, lddrepair,subcatchment, report, pcr2numpy, boolean, nominal
import numpy as np

from ..global_modules.add1 import loadmap, compressArray, decompress, generateName, loadLAI
from ..global_modules.settings import MaskInfo, LisSettings
from . import HydroModule


class evapowater(HydroModule):
    """
    # ************************************************************
    # ***** EVAPORATION FROM OPEN WATER **************************
    # ************************************************************
    """
    input_files_keys = {
        'openwaterevapo': ['LakeMask', 'maxNoEva'],
        'varfractionwater': ['FracMaxWater', 'WFractionMaps']
    }
    module_name = 'EvapoWater'

    def __init__(self, evapowater_variable):
        self.var = evapowater_variable

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------

    def initial(self):
        """ initial part of the evapo water module
        """

        # ************************************************************
        # ***** EVAPORATION
        # ************************************************************
        self.var.EvaCumM3 = MaskInfo.instance().in_zero()
        self.var.EvaWBM3 = MaskInfo.instance().in_zero()
        # water use cumulated amount
        # water use substep amount
        settings = LisSettings.instance()
        option = settings.options
        binding = settings.binding
        maskinfo = MaskInfo.instance()
        if option['openwaterevapo']:

            # exclude all waterbodies from calculation of evaporation, because the are calculated in the reservoir, lake, wetland routine
            # to avoid double counting
            # summing up all lakes, reservoir, wetlands points
            waterbody = maskinfo.in_zero()
            # Reservoir is taken out because there are no metadata on area for reservoirs
            #if option['simulateReservoirs']:
            #    waterbody += self.var.ReservoirSitesC

            # optional use evaporation directly from waterbodies by using the area of a waterbody
            # otherwise it using the previous method (or if InitLisflood is used)
            if not(option['InitLisflood']):
                if not('openwatereva_area' in option):
                    option['openwatereva_area'] = False
                if option['simulateLakes'] and option['openwatereva_area']:
                    waterbody += self.var.LakeSitesC
                if option['simulateReservoirs'] and option['openwatereva_area']:
                    waterbody += self.var.ReservoirSitesC
                # if you use wetlands than you have to use openwatereva, because it is using variable area per day
                if option['simulateWetlands']:
                    waterbody += self.var.WetlandSitesC
            waterbody[waterbody>0] = 1
            waterbodyPcr = boolean(decompress(waterbody))
            # creating a subcatchment upstream of all waterbody points
            sub = subcatchment(self.var.LddStructuresKinematic,waterbodyPcr)
            LakeMask = loadmap('LakeMask', pcr=True)
            # exclude all waterbody with are waterboy point and lakemask cells upstream of waterbody points
            # (those lakemask cells with are not unpstream, are kept for calculation
            evaMask = (boolean(sub) & boolean(LakeMask)) | boolean(waterbodyPcr)
            # invert and compress  => every cell who should be used has a 1 , all other = 0
            # the mask is used to avoid double accounting of waterbodies for open water evapo
            self.var.evaMask = -(compressArray(evaMask)-1)

            #Using the previous method - substracting evaporation from fraction of water in a gridcell
            lmask = ifthenelse(LakeMask != 0, self.var.LddStructuresKinematic, 5)
            LddEva = lddrepair(lmask)
            lddC = compressArray(LddEva)
            inAr = decompress(np.arange(maskinfo.info.mapC[0], dtype="int32"))
            self.var.downEva = (compressArray(downstream(LddEva, inAr))).astype("int32")
            # each upstream pixel gets the id of the downstream pixel
            self.var.downEva[lddC == 5] = maskinfo.info.mapC[0]
            self.var.maxNoEva = int(loadmap('maxNoEva'))
            # all pits gets a high number


    def dynamic(self):
        """ dynamic part of the evaporation from open water
        """
        settings = LisSettings.instance()
        option = settings.options
        if option['openwaterevapo']:
            # ***********************************************
            # *********  EVAPORATION FROM OPEN WATER  *******
            # ***********************************************
            if  option['openwatereva_area'] and not option['InitLisflood']:
                # for waterbodies  tjhe evaporation is substracted pot water evaporation x waterbody area
                # for wetlands this option has to be used if you want to use changing area
                # This calculates the evaporation from open water as a fraction of each gridcell
                # but excludes waterbodies, because the are treated separately in lakes, reservoir, wetland modules.
                UpstreamEva = self.var.EWRef * self.var.MMtoM3 * self.var.WaterFraction
                # only for those where there is no reservoir, lake, wetland cell
                # evaMask is taking out the cells with are upstream of waterbodies and are in
                UpstreamEva = UpstreamEva * self.var.evaMask
                # evaporation for loop is amount of water per timestep [cu m]
                # Volume of potential evaporation from water surface  per time step (conversion to [m3])
                ChanMIter = self.var.ChanM3Kin.copy()
                # for Iteration loop: First value is amount of water in the channel
                # amount of water in bankful (first line of routing)
                self.var.EvaAddM3 = np.minimum(UpstreamEva, ChanMIter * 0.9)
                # 10% of the discharge must stay in the river


            else:
                # This calculates the evaporation from open water as a fraction of each gridcell,
                # also for waterbodies
                # For waterbodies it is trying to collect water from upstream cells (default of maxNoEva=5)
                UpstreamEva = self.var.EWRef * self.var.MMtoM3 * self.var.WaterFraction
                # evaporation for loop is amount of water per timestep [cu m]
                # Volume of potential evaporation from water surface  per time step (conversion to [m3])
                ChanMIter = self.var.ChanM3Kin.copy()
                # for Iteration loop: First value is amount of water in the channel
                # amount of water in bankful (first line of routing)
                ChanLeft = ChanMIter * 0.1
                # 10% of the discharge must stay in the river
                self.var.EvaAddM3 = MaskInfo.instance().in_zero()
                #   real water consumption is set to 0

                for NoEvaExe in range(self.var.maxNoEva):
                    ChanHelp = np.maximum(ChanMIter - UpstreamEva, ChanLeft)
                    EvaIter = np.maximum(UpstreamEva - (ChanMIter - ChanHelp), 0)
                    # new amount is amout - evaporation use till a limit
                    # new evaporation is evaporation - water is used from channel network
                    ChanMIter = ChanHelp.copy()
                    self.var.EvaAddM3 += UpstreamEva - EvaIter
                    # evaporation is added up; the sum is the same as sum of original water use
                    # UpstreamEva = upstream(self.var.LddEva,EvaIter)
                    UpstreamEva = np.bincount(self.var.downEva, weights=EvaIter)[:-1]
                    # remaining water use is moved down the the river system,

            self.var.EvaAddM3Dt = self.var.EvaAddM3 * self.var.InvNoRoutSteps
            # splitting water use per timestep into water use per sub time step
            self.var.EvaCumM3 += self.var.EvaAddM3
            self.var.EvaWBM3 = self.var.EvaAddM3
            # summing up for water balance calculation
