

# keys = ['filename', 'name', 'smiles', 'abs_max', 'abs_lambda_max', 'abs_lambda_end', 'uv_reference', 'uv_absorption', 'uv_maintenance', \
#         'PL_max', 'PL_lambda_max', 'PL_maintenance', 'PL_rgb', 'relative_QY', 'max_gain_factor', 'max_gain_wavelength', 'PL_lifetime', \
#         'radiative_decay', 'Gain_cross_section', 'Gain*PL_retention']

config_dict = {
    'Gain_cross_section' : {'Label' : 'Gain_CS', 'Normalize_at' : 4E-16, 'unit' : 'cm2'},
    'relative_QY' :  {'Label' : 'PLQY', 'Normalize_at' : 1, 'unit' : ''},
    'radiative_decay' :  {'Label' : 'radiative_decay', 'Normalize_at' : 1.2E9, 'unit' : 'ns-1'}, 
    'max_gain_factor' :  {'Label' : 'Spectral_factor', 'Normalize_at' : 2E-24, 'unit' : 'cm2 s'},
    'PL_maintenance' :  {'Label' : 'PL_retention', 'Normalize_at' : 1, 'unit' : ''},
    'Gain*PL_retention' : {'Label' : 'Gain_CS*PL_retention', 'Normalize_at' : 4E-16, 'unit' : 'cm2'}
    
}

# d = [krd * Spectral_factor * PL_retention, PLQY, krd, Spectral_factor, PL_retention,  krd * Spectral_factor]