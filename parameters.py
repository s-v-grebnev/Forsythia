"""
(c) De Feo, Tulebaev, Costello, Grebnev
"""

parameters = {
    'gost256': {'lA': 2, 'lB': 3, 'eA': 372, 'eB': 239, 'f': 7, 'pm1': -1},
    'gost128_1': {'lA': 2, 'lB': 3, 'eA': 208, 'eB': 129, 'f': 5, 'pm1': -1},
    'gost128_2': {'lA': 2, 'lB': 3, 'eA': 253, 'eB': 161, 'f': 7, 'pm1': -1},
    'limonnitsa': {'lA': 2, 'lB': 3, 'eA': 451, 'eB': 284, 'f': 1, 'pm1': -1},
    'sike434': {'lA': 2, 'lB': 3, 'eA': 216, 'eB': 137, 'f': 1, 'pm1': -1},
    'sike503': {'lA': 2, 'lB': 3, 'eA': 250, 'eB': 159, 'f': 1, 'pm1': -1},
    'sike610': {'lA': 2, 'lB': 3, 'eA': 305, 'eB': 192, 'f': 1, 'pm1': -1},
    'sike751': {'lA': 2, 'lB': 3, 'eA': 372, 'eB': 239, 'f': 1, 'pm1': -1},
    '2-3-8': {'lA': 2, 'lB': 3, 'eA': 6, 'eB': 1, 'f': 1, 'pm1': -1},
    '2-3-40': {'lA': 2, 'lB': 3, 'eA': 22, 'eB': 15, 'f': 1, 'pm1': -1},
}

# Alice basis (304*z + 361 : 414*z + 157 : 113*z + 299), (368*z + 34 : 124*z + 103 : 212*z + 203)
# Bob's basis (172*z + 176 : 319*z + 106 : 70*z), (154*z + 26 : 18*z + 281 : 298*z + 416))

toy_example = {'lA': 2, 'lB': 3, 'eA':  4, 'eB': 3, 'f': 1, 'A': [104, 0], 'C': [1,0],
               'xp2': [361, 304], 'yp2': [157, 414], 'zp2': [299, 113],
               'xq2': [34, 368], 'yq2': [103, 124], 'zq2': [203, 212],
               'xp3': [176, 172], 'yp3': [106, 319], 'zp3': [0, 70],
               'xq3': [26, 154], 'yq3': [281, 18], 'zq3': [416, 298]
               }

costello = {'lA': 2, 'lB': 3, 'eA':  4, 'eB': 3, 'f': 1, 'A': [423, 329], 'C': [1,0],
            'xp2': [248, 100], 'yp2': [199, 304], 'zp2': [1, 0],
            'xq2': [394, 426], 'yq2': [79, 51], 'zq2': [1, 0],
            'xp3': [275, 358], 'yp3': [104, 410], 'zp3': [1, 0],
            'xq3': [185, 20], 'yq3': [239, 281], 'zq3': [1,0]
            }