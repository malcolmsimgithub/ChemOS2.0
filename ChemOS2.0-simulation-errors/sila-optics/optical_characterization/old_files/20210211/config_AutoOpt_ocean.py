


refractive_indices = {
	'toluene' : 1.5062,
	'ACN' : 1.3467,
	'Ethanol' : 1.3657
}


config = {
	"absorption": {
		"abs_minimum_volume" : 0.09,
		"abs_maximum_absorption" : 1,
        "abs_equibliration_time" :30,
        "abs_exposure" : 0.14,
		"dark_average" : 30,
		"abs_average" : 120,
		"abs_draw_velocity" : 50,
		"abs_dispense_velocity" : 300,
		"abs_do_plot" : False,
		"abs_calc_range" : [300,800]
	},

	"PL": {
		"PL_minimum_volume" : 0.12,
		"PL_equibliration_time" : 30,
		"led_power" : 30,
        "PL_average" : 30,
        "PL_max_exposure" : 1,
		"PL_min_measurement_time" : 10,
        "PL_initial_exposure" : 0.3,
        "PL_target_intensity" : 0.4,
        "excitation_wavelength" : "365 nm",
        "uv_average" : 10,
		"PL_draw_velocity" : 50, 
		"PL_dispense_velocity" : 300,
		"PL_calc_range" : [380,900],
		"PL_cutoff" : 0.95
	},

	"TE" : {
		"TE_minimum_volume" : 0.12,
		"TE_draw_velocity" : 500,
		"TE_dispense_velocity" : 1000,
	    "min_rate" : 1E-2, 
        "max_rate" : 3E-2, 
        "accumulation" : 10000,
        "initial_frequency" : 1E6,
		"initial_filter_position" : 12,
		"fit_order" : 1,
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
		"N_parallel" : 1,
		"waiting_timeout" : 600
	},

	"data_processing": {
        "filter_size": 5,
        "correction_datafile" : "C:/Users/MatterLab/Dropbox/PythonScript/kazu/optical_characterizations/20210201_spectral_response_curve.txt",
		"absorption_threshold" : -0.0005,
		"quantum_yeild_reference" : {'value' : 2.6E5, 'QY' : 0.89, 'solvent' : 'toluene'}   #https://pubs.acs.org/doi/10.1021/acsmaterialslett.9b00536
    }, 

	"data_path": {
        "job_input": "C:/Users/MatterLab/Dropbox/PythonScript/HPLCMS_characterization/sample_to_measure",
		"job_processing": "C:/Users/MatterLab/Dropbox/PythonScript/HPLCMS_characterization/sample_waiting/",
        "status" : "C:/Users/MatterLab/Dropbox/PythonScript/HPLCMS_characterization/status"
    } 
}
