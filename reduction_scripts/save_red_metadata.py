import pickle

#------------------------------------------------------
bias_obs = ['r0602','r0607','r0610','r0612','r0615']

domeflat_obs = ['r0603']

twiflat_obs = ['r0412','r0413','r0414']

dark_obs = []

arc_obs = ['r0608']

wire_obs = ['r0605','r0313','r0314']

#------------------------------------------------------
sci_obs = [
    # HCG91-X (625)
    #{'sci'  : ['r0625'],
    # 'arc'  : ['r0621'],
    # 'sky'  : ['r0626'],
    #'bias' : ['r0623']
    # },
    # HCG91-X (627)
    #{'sci'  : ['r0627'],
    # 'arc'  : ['r0630'],
    # 'sky'  : ['r0626'],
    #'bias' : ['r0632']
    # },
    # HCG91-X (633)
    #{'sci'  : ['r0633'],
    # 'arc'  : ['r0630'],
    # 'sky'  : ['r0634'],
    #'bias' : ['r0632']
    # },
    # HCG91-X (635)
    #{'sci'  : ['r0635'],
    # 'arc'  : ['r0637'],
    # 'sky'  : ['r0634'],
    #'bias' : ['r0639']
    # },
    # HCG91-D (640)
    #{'sci'  : ['r0640'],
    # 'arc'  : ['r0637'],
    # 'sky'  : ['r0641'],
    # 'bias' : ['r0639']
    # },
    # HCG91-D (642)
    #{'sci'  : ['r0642'],
    # 'arc'  : ['r0645'],
    # 'sky'  : ['r0641'],
    # 'bias' : ['r0647']
    # },
    # HCG91-D (648)
    {'sci'  : ['r0648'],
     'arc'  : ['r0645'],
     'sky'  : ['r0649'],
     'bias' : ['r0647']
     },
    # HCG91-D (650)
    #{'sci'  : ['r0650'],
    # 'arc'  : ['r0652'],
    # 'sky'  : ['r0649'],
    # 'bias' : ['r0655']
    # },


    ]

#------------------------------------------------------
std_obs = [   
    # EG131 (624)
    #{'sci'  : ['r0624'],
    # 'sky'  : [],
    # 'arc'  : ['r0621'],
    # 'bias' : ['r0623'],
    # 'type' : ['telluric']
    # },
    # HD204543 (629)
    #{'sci'  : ['r0629'],
    # 'sky'  : [],
    # 'arc'  : ['r0630'],
    # 'bias' : ['r0632'],
    # 'type' : ['flux', 'telluric']
    # },
    # Feige110 (636)
    {'sci'  : ['r0636'],
     'sky'  : [],
     'arc'  : ['r0637'],
     'bias' : ['r0639'],
     'type' : ['flux', 'telluric']
     },
    # Feige110 (643)
    #{'sci'  : ['r0643'],
    # 'sky'  : [],
    # 'arc'  : ['r0645'],
    # 'bias' : ['r0647'],
    # 'type' : ['flux', 'telluric']
    # },
    # Feige110 (643)
    #{'sci'  : ['r0651'],
    # 'sky'  : [],
    # 'arc'  : ['r0652'],
    # 'bias' : ['r0655'],
    # 'type' : ['flux', 'telluric']
    # },
    

]

#------------------------------------------------------
night_data = {
    'bias' : bias_obs,
    'domeflat' : domeflat_obs,
    'twiflat' : twiflat_obs,
    'dark' : dark_obs,
    'wire' : wire_obs,
    'arc'  : arc_obs,
    'sci'  : sci_obs,
    'std'  : std_obs}

f1 = open('wifesR_20120819_metadata.pkl', 'w')
pickle.dump(night_data, f1)
f1.close()
