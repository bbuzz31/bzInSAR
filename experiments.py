EXP_DICT  = dict(root='HR',                     # root path (location)
                  naz=0, nrg=0,                 # number of azimuth, range looks
                  nx=0,  ny=0,                  # window size
                  sref     ='SR',               # use the single reference stack
                  ann      =False,              # use annual pairs
                  sbas     =False,              # use small baseline stack
                  gacos    =False,              # apply GACOS corrections
                  era5     =False,              # apply GACOS corrections
                  gap_fill =True,               # gap fill prior to unwrapping
                  bridging =False,              # apply bridging error correct
                  closure  =False,              # apply phase closure correction
                  cust_net ='',                 # custom network (specify ifgs)
                  custom   ='',                 # custom ifg ext, e.g., _sum
                  satellite='Sentinel1'
                  )


## for making IFGS with stack processor

HR_SR  = {**EXP_DICT,  'nx': 33, 'ny': 15, 'naz':2, 'nrg':5} # until 2023
# HR_SR0 = {**HR_SR, 'root':'HR0'} # until 2022; deleted
HR_SR2 = {**HR_SR,  'nx': 55, 'ny': 25}

# mostly just for unwrapping
HR_SBAS  = {**HR_SR, 'sref': False, 'sbas': 3}
HR_ANN   = {**HR_SR, 'ann': 2}

HR_2016  = {**HR_SR, 'root': 'HR_2016', 'nx': 11, 'ny': 5}

HR_Gacos = {**HR_SR2, 'gacos':True}
HR_Era5  = {**HR_SR2, 'era5':True}

HR_ANNW     = {**HR_SR2, 'ann': '12-1-2'} # months
HR_ANNW_90  = {**HR_ANNW, 'sref': False, 'naz': 7, 'nrg': 19}

HR_SBAS1    = {**HR_SR2, 'sref': 'SR', 'sbas': 2}  # tests
HR_SBAS2    = {**HR_SR2, 'sref': 'SR', 'ann': '12-1-2', 'sbas': 2}   # tests
HR_SBAS3    = {**HR_SR2, 'sref': False, 'ann': '12-1-2', 'sbas': 2}   # tests
HR_SBAS_90  = {**HR_SBAS, 'naz': 7, 'nrg': 19}
HR_SBAS1_90  = {**HR_SBAS1, 'naz': 7, 'nrg': 19}

# unwrapping corrections need some redundancy
HR_ANNUP  = {**HR_ANN, 'closure':True}
HR_ANNUB  = {**HR_ANN, 'bridging':True} # need to mkae timeseries
HR_ANNUPB = {**HR_ANN, 'closure':True, 'bridging':True} # to make in stack


HR_FR = {**HR_SR2, 'naz':1, 'nrg':1, 'gap_fill':False, 'custom':'_sum'}
HR_90 = {**HR_SR2, 'naz':7, 'nrg':19}

## custom pairs
HR_90c = {**HR_SR2, 'naz':7, 'nrg':19, 'sref': False, 'cust_net': '2020'}

HR_90i = {**HR_90, 'sref':'SRi'} # invert single ref network
HR3_90 = {**HR_90, 'root': 'HR3'}

HR_ARIA = {**EXP_DICT, 'root': 'HR_ARIA', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2}
HR_ARIA_24 = {**EXP_DICT, 'root': 'HR_ARIA_2024', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2} # new products, new AT
HR_ARIA_24a = {**EXP_DICT, 'root': 'HR_ARIA_2024a', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2} # old products, new AT

HR_SBAS20 = {**HR_90, 'sref': '', 'cust_net': '2020'}   # match 2020, ML
HR_SBAS20a = {**HR_SR2, 'sref': '', 'cust_net': '2020a'}   #
HR_SBASa   = {**HR_90, 'sref': '', 'cust_net': 'ARIA'}   # all poss ifgs

# tests
HR_Test = {**HR_FR, 'root':'HR_Test3'}


# -------------------------------- NYC---------------------------------------- #
NYC_SR0  = {**EXP_DICT, 'root':'NYC', 'nx':11, 'ny':5, 'naz':2, 'nrg':5}
NYC_SRw  = {**EXP_DICT, 'root':'NYC', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
NYC_SR   = {**EXP_DICT, 'root':'NYC1', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
NYC_SRc  = {**EXP_DICT, 'root':'NYC2', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
NYC_SRd  = {**EXP_DICT, 'root':'NYC_Dolphin', 'nx':33, 'ny':15, 'naz':1, 'nrg':1}

NYC_SBAS = {**NYC_SR, 'sref': False, 'sbas': 3}
NYC_ANN  = {**NYC_SR, 'ann': 2}
NYC_ARIA  = {**EXP_DICT, 'root':'NYC_ARIA', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2}
NYC_ARIA1 = {**EXP_DICT, 'root':'NYC_ARIA1', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2} # bigger region

NYC_ARIA_ATM0 = {**EXP_DICT, 'root':'NYC_ARIA_ATM', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2} # slant tropo corrs
NYC_ARIA_ATM  = {**EXP_DICT, 'root':'NYC_ARIA_HRES', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2} # slant tropo corrs

NYC_ARIA_SRTM = {**EXP_DICT, 'root':'NYC_ARIA_SRTM', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2}

NYC_ALOS    = {**EXP_DICT, 'root':'NYC', 'naz':9, 'nrg':4, 'sbas':3, 'satellite': 'ALOS'}
NYC_ALOS_FR = {**EXP_DICT, 'root':'NYC_FRInGE', 'naz':1, 'nrg':1, 'satellite': 'ALOS'}


# -------------------------------- NJ ---------------------------------------- #
NNJ_SR = {**EXP_DICT, 'root': 'NNJ', 'nx':33, 'ny': 15, 'naz':2, 'nrg':5}

NJ_SR = {**EXP_DICT, 'root': 'NJ', 'nx':33, 'ny': 15, 'naz':2, 'nrg':5}


# -------------------------------- SC---------------------------------------- #
# SC_Base = {**EXP_DICT, 'root':'SC', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}

# SC_SBAS2 = {**SC_Base, 'sref': 'SR', 'ann': '12-1-2', 'sbas': 2}   # tests

# SC_SBASa = {**SC_Base, 'sref': '', 'cust_net': 'ARIA'}   # tests

# Charleston_SR   = {**EXP_DICT, 'root':'Charleston', 'nx':33, 'ny':15, 'naz':2, 'nrg':5} ## deleted
Charleston_SR2  = {**EXP_DICT, 'root':'Charleston2', 'nx':33, 'ny':15, 'naz':2, 'nrg':5} # multiple frames to 2023
Charleston_SR_Small = {**EXP_DICT, 'root': 'Charleston_Small1', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
Charleston_SR_Test  = {**EXP_DICT, 'root':'Charleston_TEST', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
Charleston_SBAS = {**Charleston_SR2, 'sref': False, 'sbas': 3}
Charleston_ANN  = {**Charleston_SR2, 'ann': 2}

Charleston_SBASa = {**Charleston_SR2, 'sref': '', 'cust_net': 'ARIA'}   # tests
# Charleston_SBAS2 = {**Charleston_Base, 'sref': 'SR', 'ann': '12-1-2', 'sbas': 2}   # tests

Charleston_ARIA = {**EXP_DICT, 'root': 'Charleston_ARIA', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2}

# Savannah_SR   = {**EXP_DICT, 'root':'Savannah', 'nx':33, 'ny':15, 'naz':2, 'nrg':5} # deleted
Savannah_SR2  = {**EXP_DICT, 'root':'Savannah2', 'nx':33, 'ny':15, 'naz':2, 'nrg':5} # multiple frames to 2023
Savannah_SBAS = {**Savannah_SR2, 'sref': False, 'sbas': 3}
Savannah_ANN  = {**Savannah_SR2, 'ann': 2}
Savannah_ARIA = {**EXP_DICT, 'root': 'Savannah_ARIA', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2}


# -------------------------------- CA---------------------------------------- #
SF_SR   = {**EXP_DICT, 'root':'SF', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
LongBeach_ARIA = {**EXP_DICT, 'root': 'LongBeach_ARIA', 'naz':7, 'nrg':19, 'sref':False, 'sbas':2}

# -------------------------------- FL---------------------------------------- #
Miami_SR   = {**EXP_DICT, 'root':'Miami', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
Miami_SRc  = {**EXP_DICT, 'root':'Miami1', 'nx':33, 'ny':15, 'naz':2, 'nrg':5} #
Kennedy_SR  = {**EXP_DICT, 'root':'Kennedy', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
Orlando_SR  = {**EXP_DICT, 'root':'Orlando', 'nx':33, 'ny':15, 'naz':2, 'nrg':5} # symlinked with Kennedy

# -------------------------------- PA --------------------------------------- #
Philly_SR  = {**EXP_DICT, 'root':'Philly', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}

# ------------------------------ Texas -------------------------------------- #
Houston_SR  = {**EXP_DICT, 'root':'Houston', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}

# ------------------------------ DC -------------------------------------- #
DC_SR  = {**EXP_DICT, 'root':'DC', 'nx':33, 'ny':15, 'naz':2, 'nrg':5}
DC_ARIA  = {**EXP_DICT, 'root':'DC_ARIA', 'nx':7, 'ny':19, 'sref': False, 'sbas': 2}


# ----------------------------- Islands ------------------------------------- #
## Stripmap;
# Kiribati SR0 is complete and good but uses synthetic ifgs; use all instaed
# Kiribati SR is complete and use TG-ALT reference
Kiribati_SR0 = {**EXP_DICT, 'root':'Kiribati0', 'nx':33, 'ny':15, 'naz':3, 'nrg':3}
Kiribati_SR  = {**EXP_DICT, 'root':'Kiribati', 'nx':33, 'ny':15, 'naz':3, 'nrg':3}
