import pickle

#------------------------------------------------------
bias_obs = [
   'OBK-239552-WiFeS-Red--UT20230306T055451-8',
   'OBK-239552-WiFeS-Red--UT20230306T055354-6',
    ]

domeflat_obs = [
   'OBK-239488-WiFeS-Red--UT20230306T054004-0',
   'OBK-239360-WiFeS-Red--UT20230306T054752-8',
   'OBK-239424-WiFeS-Red--UT20230306T054636-7',
   'OBK-239360-WiFeS-Red--UT20230306T055257-3',
   'OBK-239488-WiFeS-Red--UT20230306T053542-0',
   'OBK-239360-WiFeS-Red--UT20230306T054908-9',
   'OBK-239520-WiFeS-Red--UT20230306T053858-3',
   'OBK-239488-WiFeS-Red--UT20230306T054109-3',
   'OBK-239520-WiFeS-Red--UT20230306T053753-0',
   'OBK-239488-WiFeS-Red--UT20230306T053647-3',
   'OBK-239456-WiFeS-Red--UT20230306T054320-3',
   'OBK-239360-WiFeS-Red--UT20230306T055141-2',
   'OBK-239424-WiFeS-Red--UT20230306T054531-3',
   'OBK-239424-WiFeS-Red--UT20230306T054426-0',
   'OBK-239360-WiFeS-Red--UT20230306T055025-1',
   'OBK-239456-WiFeS-Red--UT20230306T054215-0',
    ]

twiflat_obs = [
   'OBK-14240-WiFeS-Red--UT20230305T083835-2',
   'OBK-14240-WiFeS-Red--UT20230304T083949-9',
   'OBK-14240-WiFeS-Red--UT20230306T083750-0',
    ]

dark_obs = [
    ]

arc_obs = [
   'OBK-239392-WiFeS-Red--UT20230306T053228-2',
   'OBK-239392-WiFeS-Red--UT20230306T053340-8',
   'OBK-239392-WiFeS-Red--UT20230306T053115-7',
    ]

wire_obs = [
    ]

#------------------------------------------------------
sci_obs = [
    # GRB230307A
    {'sci'  : ['OBK-245120-WiFeS-Red--UT20230310T105634-4',
               'OBK-245120-WiFeS-Red--UT20230310T101435-8',
               'OBK-245120-WiFeS-Red--UT20230310T103535-1'],
     'sky'  : []},
    ]

#------------------------------------------------------
std_obs = [
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

f1 = open('wifesR_20230312_metadata.pkl', 'wb')
pickle.dump(night_data, f1)
f1.close()
