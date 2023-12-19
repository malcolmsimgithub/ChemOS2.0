#!/usr/bin/env 

#========================================================================

import sys
import json
import numpy as np 
#========================================================================

class Configuration(object):

	def __init__(self, me = ''):
		self.me          = me + ':'
		self.added_props = []
		self.added_attrs = []

	def __str__(self):
		new_line = '%s\n' % self.me
		for prop in sorted(self.added_props):
			new_line += '--> %s:\t%s\n' % (prop, getattr(self, prop))
		return(new_line)


	def __iter__(self):
		for _ in range(self.num_elements):
			info_dict = {}
			for prop_index, prop in enumerate(self.added_props):
				info_dict[prop] = self.added_attrs[prop_index][_]
			yield info_dict


	def __getitem__(self, index):
		info_dict = {}
		for prop_index, prop in enumerate(self.added_props):
			info_dict[prop] = self.added_attrs[prop_index][index]
		return info_dict

	
	def to_dict(self):
		return {prop: getattr(self, prop) for prop in sorted(self.added_props)}


	def add_attr(self, prop, attr):
		setattr(self, prop, attr)
		if not prop in self.added_props:		
			self.added_props.append(prop)
			self.added_attrs.append(attr)
		try:
			self.num_elements = len(attr)
		except TypeError:
			pass

	def get_attr(self, prop):
		return getattr(self, prop)

#========================================================================

class ConfigParser():

	def __init__(self, default, **kwargs):

		self.config_dict = default
		self.config_file = None

		for key1, val1 in kwargs.items():
			for key2, val2 in val1.items():
				self.config_dict[key1][key2] = val2

	def _parse(self, provided_settings):

		self.data_path = Configuration('data_path')
		for key, value in provided_settings['data_path'].items():
			self.data_path.add_attr(key, value)

		self.absorption = Configuration('absorption')
		for key, value in provided_settings['absorption'].items():
			self.absorption.add_attr(key, value)

		self.PL = Configuration('PL')
		for key, value in provided_settings['PL'].items():
			self.PL.add_attr(key, value)

		self.TE = Configuration('TE')
		for key, value in provided_settings['TE'].items():
			self.TE.add_attr(key, value)

		self.dilution = Configuration('dilution')
		for key, value in provided_settings['dilution'].items():
			self.dilution.add_attr(key, value)

		self.redissolution = Configuration('redissolution')
		for key, value in provided_settings['redissolution'].items():
			self.redissolution.add_attr(key, value)

		self.data_processing = Configuration('data_processing')
		for key, value in provided_settings['data_processing'].items():
			self.data_processing.add_attr(key, value)

		#self.peak_search = Configuration('peak_search')
		# self.peak_search.add_attr('DAD', provided_settings['DAD_peak_search_setting'])
		# self.peak_search.add_attr('MS_TIC', provided_settings['MS_peak_search_setting_TIC'])
		# self.peak_search.add_attr('MS_XIC', provided_settings['MS_peak_search_setting_XIC'])



	@property
	def abs_minimum_volume(self):
		return self.absorption.abs_minimum_volume

	@property
	def abs_maximum_absorption(self):
		return self.absorption.abs_maximum_absorption

	@property
	def abs_equibliration_time(self):
		return self.absorption.abs_equibliration_time

	@property
	def abs_exposure(self):
		return self.absorption.abs_exposure

	@property
	def abs_dark_average(self):
		return self.absorption.dark_average

	@property
	def abs_average(self):
		return self.absorption.abs_average

	@property
	def abs_draw_velocity(self):
		return self.absorption.abs_draw_velocity

	@property
	def abs_dispense_velocity(self):
		return self.absorption.abs_dispense_velocity

	@property
	def abs_do_plot(self):
		return self.absorption.abs_do_plot

	@property
	def abs_calc_range(self):
		return self.absorption.abs_calc_range

	@property
	def PL_minimum_volume(self):
		return self.PL.PL_minimum_volume

	@property
	def PL_equibliration_time(self):
		return self.PL.PL_equibliration_time

	@property
	def PL_led_power(self):
		return self.PL.led_power

	@property
	def PL_average(self):
		return self.PL.PL_average

	@property
	def PL_max_exposure(self):
		return self.PL.PL_max_exposure

	@property
	def PL_min_measurement_time(self):
		return self.PL.PL_min_measurement_time

	@property
	def PL_initial_exposure(self):
		return self.PL.PL_initial_exposure

	@property
	def PL_target_intensity(self):
		return self.PL.PL_target_intensity

	@property
	def PL_excitation_wavelength(self):
		return self.PL.excitation_wavelength

	@property
	def PL_uv_average(self):
		return self.PL.uv_average

	@property
	def PL_draw_velocity(self):
		return self.PL.PL_draw_velocity

	@property
	def PL_dispense_velocity(self):
		return self.PL.PL_dispense_velocity

	@property
	def PL_calc_range(self):
		return self.PL.PL_calc_range

	@property
	def PL_cutoff(self):
		return self.PL.PL_cutoff

	@property
	def TE_minimum_volume(self):
		return self.TE.TE_minimum_volume

	@property
	def TE_draw_velocity(self):
		return self.TE.TE_draw_velocity

	@property
	def TE_dispense_velocity(self):
		return self.TE.TE_dispense_velocity

	@property
	def TE_min_rate(self):
		return self.TE.min_rate

	@property
	def TE_max_rate(self):
		return self.TE.max_rate

	@property
	def TE_accumulation(self):
		return self.TE.accumulation
		
	@property
	def TE_initial_frequency(self):
		return self.TE.initial_frequency

	@property
	def TE_initial_filter_position(self):
		return self.TE.initial_filter_position

	@property
	def TE_fit_order(self):
		return self.TE.fit_order

	@property
	def TE_fit_weight(self):
		return self.TE.fit_weight

	@property
	def dilution_draw_velocity(self):
		return self.dilution.dilution_draw_velocity

	@property
	def dilution_dispense_velocity(self):
		return self.dilution.dilution_dispense_velocity

	@property
	def evap_temp(self):
		return self.redissolution.evaporation_temparature

	@property
	def evap_time(self):
		return self.redissolution.evaporation_time

	@property
	def dissolution_temp(self):
		return self.redissolution.dissolution_temparature

	@property
	def dissolution_time(self):
		return self.redissolution.dissolution_time

	@property
	def shaking_speed(self):
		return self.redissolution.shaking_speed

	@property
	def dissolution_volume(self):
		return self.redissolution.solvent_volume

	@property
	def dissolution_idling_temp(self):
		return self.redissolution.idling_temparature

	@property
	def dissolution_waiting_timeout(self):
		return self.redissolution.waiting_timeout

	@property
	def filter_size(self):
		return self.data_processing.filter_size

	@property
	def correction_datafile(self):
		return self.data_processing.correction_datafile
	
	@property
	def absorption_threshold(self):
		return self.data_processing.absorption_threshold

	@property
	def quantum_yeild_reference(self):
		return self.data_processing.quantum_yeild_reference

	@property
	def data_path_job_input(self):
		return self.data_path.job_input

	@property
	def data_path_job_processing(self):
		return self.data_path.job_processing

	@property
	def data_path_status(self):
		return self.data_path.status

	# @safe_execute(GryffinParseError)
	def parse_config_file(self, config_file = None):

		if not config_file is None:
			self.config_file = config_file

		self._parse(self.config_dict)


	# @safe_execute(GryffinParseError)
	def parse_config_dict(self, config_dict = None):
	
		if not config_dict is None:
			self.config_dict = config_dict

		self._parse(self.config_dict)


	def parse(self):
		# test if both dict and file have been provided
		if self.config_dict is not None and self.config_file is not None:
			print('Found both configuration file and configuration dictionary. Will parse configuration from dictionary and ignore file', 'WARNING')
			self.parse_config_dict(self.config_dict)
		elif self.config_dict is not None:
			self.parse_config_dict(self.config_dict)
		elif self.config_file is not None:
			raise NotImplementedError('Loading settings from a config file is not implemented yet. Please use config dict.')
		else:
			print('Cannot parse configuration due to missing configuration file or configuration dictionary', 'ERROR')



	
