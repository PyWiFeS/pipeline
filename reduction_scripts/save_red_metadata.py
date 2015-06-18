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

    ]

#------------------------------------------------------
std_obs = [   
    # Feige110 (636)
    {'sci'  : ['r0636'],
     'sky'  : [],
     'name' : [],
     'arc'  : ['r0637'],
     'bias' : ['r0639'],
     'type' : ['flux', 'telluric']
     },

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
