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
    # HCG91-X (627)
    #{'sci'  : ['b0627'],
    # 'arc'  : ['b0631'],
    # 'sky'  : ['b0626'],
    # 'bias' : ['b0632']
    # },
    # HCG91-X (633)
    #{'sci'  : ['b0633'],
    # 'arc'  : ['b0631'],
    # 'sky'  : ['b0634'],
    # 'bias' : ['b0632']
    # },
    # HCG91-X (635)
    #{'sci'  : ['b0635'],
    # 'arc'  : ['b0638'],
    # 'sky'  : ['b0634'],
    # 'bias' : ['b0639']
    # },
    # HCG91-D (640)
    #{'sci' : ['b0640'],
    # 'arc':['b0638'],
    # 'sky':['b0641'],
    # 'bias':['b0639']
    # }, 
    # HCG91-D (642)
    #{'sci'  : ['b0642'],
    # 'arc'  : ['b0646'],
    # 'sky'  : ['b0641'],
    # 'bias' : ['b0647']
    # },
    # HCG91-D (648)
    #{'sci'  : ['b0648'],
    # 'arc'  : ['b0646'],
    # 'sky'  : ['b0649'],
    # 'bias' : ['b0647']
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
    # HD204543 (629)
    #{'sci'  : ['b0629'],
    # 'sky'  : [],
    # 'arc'  : ['b0630'],
    # 'bias' : ['b0632'],
    # 'type' : ['flux']
    # },
    # F110 (636)
    #{'sci'  : ['b0636'],
    # 'arc'  : ['b0638'],
    # 'sky' : [],
    # 'bias': ['b0639'],
    # 'type': ['flux']
    # },
    # F110 (643)
    #{'sci'  : ['b0643'],
    # 'arc'  : ['b0646'],
    # 'sky' : [],
    # 'bias': ['b0647'],
    # 'type': ['flux']
    # },
    # Feige110 (651)
    {'sci'  : ['b0651'],
     'sky'  : [],
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
