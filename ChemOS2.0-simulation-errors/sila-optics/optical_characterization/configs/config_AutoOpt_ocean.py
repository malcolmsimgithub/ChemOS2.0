path = r"C:/Users/MatterLab/Dropbox (Aspuru-Guzik Lab)/PythonScript"
ml = 1e-3


###the refractive indices of the solvents averaged over the visible spectral region
refractive_indices = {
	'toluene' : 1.5062,
	'ACN' : 1.3467,
	'Ethanol' : 1.3657,
	'THF' : 1.4072
}


config = {
	"absorption": {
		"abs_minimum_volume" : 0.09, #in mL
		"abs_maximum_absorption" : 1,
        "abs_equilibration_time" :15,
        "abs_exposure" : 0.14,
		"dark_average" : 30,
		"abs_average" : 120,
		"abs_draw_velocity" : 300,
		"abs_dispense_velocity" : 300,
		"abs_do_plot" : False,
		"abs_calc_range" : [300, 750]
	},

	"PL": {
		"PL_minimum_volume" : 0.11,   #in mL
		"PL_equilibration_time" : 15,
		"led_power" : 30,
        "PL_average" : 30,
        "PL_max_exposure" : 1,
		"PL_min_measurement_time" : 10,
        "PL_initial_exposure" : 0.3,
        "PL_target_intensity" : 0.4,
        # "excitation_wavelength" : "365 nm",
        "uv_average" : 10,
		"PL_draw_velocity" : 300, 
		"PL_dispense_velocity" : 300,
		"PL_calc_range" : [380, 750],
		"PL_cutoff" : 0.95
	},

	"TE" : {
		"TE_minimum_volume" : 0.11,   #in mL
		"TE_draw_velocity" : 500,
		"TE_dispense_velocity" : 1000,
	    "min_rate" : 5E-3, 
        "max_rate" : 3E-2, 
        "accumulation" : 10000,
        "initial_frequency" : 1E6,
		"initial_filter_position" : 12,
		"fit_order" : 'auto',
		"fit_weight" : 0.7
	},

	"dilution": {
		"dilution_draw_velocity" : 500,
		"dilution_dispense_velocity" : 500
	},

	"redissolution":{
		"evaporation_temparature" : 50,
		"evaporation_time" : 300,
		"dissolution_temparature" : 50,
		"shaking_speed" : 8,
		"dissolution_time" : 600,
		"solvent_volume" : 0.15,
		"idling_temparature" : 50,
		"N_parallel" : 1,         #wait until this number of samples are collected in the collecton vials (only for auto_measurement with HPLCMS).
		"waiting_timeout" : 1800
	},
	
	"washing" :{
		"flow_cells_volume" : 0.5*ml, #solvent volume used for washing (per repeat) in L
		"flow_cells_repeat" : 3,
		"dilution_vial_volume" : 0.3*ml, #solvent volume used for washing (per repeat) in L
		"dilution_vial_repeat" : 3,
		"collection_vial_volume" : 0.5*ml, #solvent volume used for washing (per repeat) in L
		"collection_vial_repeat" : 3,
	},

	"data_processing": {
        "filter_size": 5,
        "correction_datafile" : "%s/kazu/optical_characterizations/calibration_curves/20210201_spectral_response_curve.txt" %path,
		"absorption_threshold" : -0.0005,
		"quantum_yeild_reference" : {'value' : 2.6E5, 'QY' : 0.89, 'solvent' : 'toluene'}   #https://pubs.acs.org/doi/10.1021/acsmaterialslett.9b00536
    }, 

	"data_path": {
        "job_input": "%s/HPLCMS_characterization/sample_to_measure" %path,
		"job_processing": "%s/HPLCMS_characterization/sample_waiting/"  %path, 
		"result_output" : "%s/HPLCMS_characterization/sample_measured/" %path,
        "status" : "%s/HPLCMS_characterization/status"  %path
    } 
}
