import pickle

#------------------------------------------------------
bias_obs = ['b0602','b0607','b0610','b0612', 'b0615']

domeflat_obs = ['b0604']

twiflat_obs = ['b0412','b0413','b0414']

dark_obs = []

arc_obs = ['b0609']

wire_obs = ['b0606', 'b0313','b0314']

#------------------------------------------------------
sci_obs = [
    # HCG91-X (625)
    #{'sci'  : ['b0625'],
    # 'arc'  : ['b0622'],
    # 'sky'  : ['b0626'],
    # 'bias' : ['b0623']
    # },
    # HCG91-D (650)
    {'sci'  : ['b0650'],
     'arc'  : ['b0653'],
     'sky'  : ['b0649'],
     'bias' : ['b0655']
     },

    ]

#------------------------------------------------------
std_obs = [    
    # Feige110 (651)
    {'sci'  : ['b0651'],
     'sky'  : [],
     #'name' : [],
     'arc'  : ['b0653'],
     'bias' : ['b0655'],
     'type' : ['flux']
     },
    
]

#------------------------------------------------------
night_data = {
    'bias' : bias_obs,
    'dark' : dark_obs,
    'domeflat' : domeflat_obs,
    'twiflat' : twiflat_obs,
    #'dark' : dark_obs,
    'wire' : wire_obs,
    'arc'  : arc_obs,
    'sci'  : sci_obs,
    'std'  : std_obs}

f1 = open('wifesB_20120819_metadata.pkl', 'w')
pickle.dump(night_data, f1)
f1.close()
