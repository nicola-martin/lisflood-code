"""

Copyright 2019-2020 European Union

Licensed under the EUPL, Version 1.2 or as soon they will be approved by the European Commission  subsequent versions of the EUPL (the "Licence");

You may not use this work except in compliance with the Licence.
You may obtain a copy of the Licence at:

https://joinup.ec.europa.eu/sites/default/files/inline-files/EUPL%20v1_2%20EN(1).txt

Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the Licence for the specific language governing permissions and limitations under the Licence.
"""

from __future__ import absolute_import
import os

import pytest


from unittest.mock import call

import lisflood
from lisflood.global_modules.add1 import loadmap
from lisflood.global_modules.errors import LisfloodError
from lisflood.main import lisfloodexe

from .test_utils import setoptions, mk_path_out


class TestOptions():
    settings_files = {'base1': os.path.join(os.path.dirname(__file__), 'data/LF_lat_lon_UseCase/run_lat_lon2.xml')}

    @classmethod
    def dummyloadmap(cls, *args, **kwargs):
        return loadmap(*args, **kwargs)


    def test_wetlands_only(self, mocker):

        settings = setoptions(self.settings_files['base1'],['simulateWetlands','reservoir_lakes_Excel','openwatereva_area'])
        mock_api = mocker.MagicMock(name='loadmap')
        mock_api.side_effect = self.dummyloadmap

        mocker.patch('lisflood.hydrological_modules.wetlands.loadmap', new=mock_api)
        lisfloodexe(settings)
        calls = [call('WetlandSites'), call('LakeMultiplier'),call('WetlandETmult'),
                call('Wetland_maxlevel'), call('WetlandInitialLevelValue'),
                call('WetlandPrevInflowValue'),call('WetlandPrevOutflowValue')]
        lisflood.hydrological_modules.wetlands.loadmap.assert_has_calls(calls)



    def test_reservoirs_only(self, mocker):
        settings = setoptions(self.settings_files['base1'],['simulateReservoirs','reservoir_lakes_Excel','openwatereva_area'])
        mock_api = mocker.MagicMock(name='loadmap')
        mock_api.side_effect = self.dummyloadmap
        mocker.patch('lisflood.hydrological_modules.reservoir.loadmap', new=mock_api)
        lisfloodexe(settings)
        calls = [call('ReservoirSites'), call('ReservoirSites', pcr=True), call('adjust_Normal_Flood'),
                 call('ReservoirRnormqMult'), call('ReservoirInitialFillValue')]
        lisflood.hydrological_modules.reservoir.loadmap.assert_has_calls(calls)
 
    def test_reservoirs_release(self, mocker):
        settings = setoptions(self.settings_files['base1'],['simulateReservoirs','reservoir_lakes_Excel',
            'reservoir_release','openwatereva_area'])
        mock_api = mocker.MagicMock(name='loadmap')
        mock_api.side_effect = self.dummyloadmap
        mocker.patch('lisflood.hydrological_modules.reservoir.loadmap', new=mock_api)
        lisfloodexe(settings)
        calls = [call('ReservoirSites'), call('ReservoirSites', pcr=True), call('adjust_Normal_Flood'),
                 call('ReservoirRnormqMult'), call('ReservoirInitialFillValue')]
        lisflood.hydrological_modules.reservoir.loadmap.assert_has_calls(calls)
        

    def test_lakes_only(self, mocker):

        settings = setoptions(self.settings_files['base1'],['simulateLakes','reservoir_lakes_Excel','openwatereva_area'])
        mock_api = mocker.MagicMock(name='loadmap')
        mock_api.side_effect = self.dummyloadmap

        mocker.patch('lisflood.hydrological_modules.lakes.loadmap', new=mock_api)
        lisfloodexe(settings)
        calls = [call('LakeSites'), call('LakeMultiplier'), call('LakeInitialLevelValue'),
                call('LakePrevInflowValue'),call('LakePrevOutflowValue')]
        lisflood.hydrological_modules.lakes.loadmap.assert_has_calls(calls)

    def test_reslakes(self, mocker):
    
        settings = setoptions(self.settings_files['base1'],['simulateLakes','simulateReservoirs','reservoir_lakes_Excel'])
        mock_api = mocker.MagicMock(name='loadmap')
        mock_api.side_effect = self.dummyloadmap

        mocker.patch('lisflood.hydrological_modules.lakes.loadmap', new=mock_api)
        lisfloodexe(settings)
        
        
    def test_reslakesWet(self, mocker):
    
        settings = setoptions(self.settings_files['base1'],['simulateLakes','simulateReservoirs','simulateWetlands',
            'reservoir_lakes_Excel','openwatereva_area',])
        mock_api = mocker.MagicMock(name='loadmap')
        mock_api.side_effect = self.dummyloadmap

        mocker.patch('lisflood.hydrological_modules.lakes.loadmap', new=mock_api)
        lisfloodexe(settings)