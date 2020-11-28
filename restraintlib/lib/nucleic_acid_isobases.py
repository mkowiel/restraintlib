

ISOCYTOSINE_PDB_CODES = ['IC']

ISOCYTOSINE_ATOM_NAMES = {
    # 'OP3': 'O3P',
    # 'P': 'P',
    # 'OP1': 'OP1',
    # 'O1P': 'OP1',
    # 'OP2': 'OP2',
    # 'O2P': 'OP2',

    # "O5'": 'O5',
    # 'O5*': 'O5',
    # "C5'": 'C5',
    # 'C5*': 'C5',
    # "C4'": 'C4',
    # 'C4*': 'C4',
    # "O4'": 'O4',
    # 'O4*': 'O4',
    # "C3'": 'C3',
    # 'C3*': 'C3',
    # "O3'": 'O3',
    # 'O3*': 'O3',
    # "C2'": 'C2',
    # 'C2*': 'C2',
    # "O2'": 'O2',
    # 'O2*': 'O2',
    # "C1'": 'C1',
    # 'C1*': 'C1',

    "N2": "N2",
    "C4": "C4",
    "N3": "N3",
    "C2": "C2",
    "O4": "O4",
    "N1": "N1",
    "C6": "C6",
    "C5": "C5",
}

ISOCYTOSINE_ATOM_RES = {
    "N2": 0,
    "C4": 0,
    "N3": 0,
    "C2": 0,
    "O4": 0,
    "N1": 0,
    "C6": 0,
    "C5": 0,
}

ISOCYTOSINE_REQUIRED_CONDITION = {
    ("N1", "C2", 2.0, 0, 0),
    ("C2", "N3", 2.0, 0, 0),
    ("N3", "C4", 2.0, 0, 0),
    ("C4", "C5", 2.0, 0, 0),
    ("C5", "C6", 2.0, 0, 0),
    ("C6", "N1", 2.0, 0, 0),
    ("C2", "N2", 2.0, 0, 0),
    ("C4", "O4", 2.0, 0, 0),
}

ISOCYTOSINE_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': [
        "aC6N1C2",
        "aN1C2N3",
        "aC2N3C4",
        "aN3C4C5",
        "aC4C5C6",
        "aC5C6N1",
        "aN2C2N1",
        "aN2C2N3",
        "aO4C4C5",
        "aO4C4N3",
    ]
}

ISOCYTOSINE_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': []
}

ISOCYTOSINE_RESTRAINTS = [
    {
        "conditions": [],
        "name": "Base==Isocytosine",
        "restraints": [
            ["dist", "dN1C2", ["N1", "C2"], 1.385, 0.009],
            ["dist", "dC2N3", ["C2", "N3"], 1.324, 0.009],
            ["dist", "dN3C4", ["N3", "C4"], 1.375, 0.009],
            ["dist", "dC4C5", ["C4", "C5"], 1.461, 0.009],
            ["dist", "dC5C6", ["C5", "C6"], 1.350, 0.009],
            ["dist", "dC6N1", ["C6", "N1"], 1.382, 0.009],
            ["dist", "dC2N2", ["C2", "N2"], 1.340, 0.009],
            ["dist", "dC4O4", ["C4", "O4"], 1.244, 0.009],

            ["angle", "aC6N1C2", ["C6", "N1", "C2"], 118.20, 0.7],
            ["angle", "aN1C2N3", ["N1", "C2", "N3"], 122.18, 0.7],
            ["angle", "aC2N3C4", ["C2", "N3", "C4"], 121.72, 0.7],
            ["angle", "aN3C4C5", ["N3", "C4", "C5"], 117.27, 0.7],
            ["angle", "aC4C5C6", ["C4", "C5", "C6"], 118.96, 0.7],
            ["angle", "aC5C6N1", ["C5", "C6", "N1"], 121.67, 0.7],
            ["angle", "aN2C2N1", ["N2", "C2", "N1"], 118.82, 0.7],
            ["angle", "aN2C2N3", ["N2", "C2", "N3"], 119.01, 0.7],
            ["angle", "aO4C4C5", ["O4", "C4", "C5"], 121.92, 0.7],
            ["angle", "aO4C4N3", ["O4", "C4", "N3"], 120.81, 0.7],
        ]
    },
]


ISOGUANINE_PDB_CODES = ['IG']

ISOGUANINE_ATOM_NAMES = {
    # 'OP3': 'O3P',
    # 'P': 'P',
    # 'OP1': 'OP1',
    # 'O1P': 'OP1',
    # 'OP2': 'OP2',
    # 'O2P': 'OP2',

    # "O5'": 'O5',
    # 'O5*': 'O5',
    # "C5'": 'C5',
    # 'C5*': 'C5',
    # "C4'": 'C4',
    # 'C4*': 'C4',
    # "O4'": 'O4',
    # 'O4*': 'O4',
    # "C3'": 'C3',
    # 'C3*': 'C3',
    # "O3'": 'O3',
    # 'O3*': 'O3',
    # "C2'": 'C2',
    # 'C2*': 'C2',
    # "O2'": 'O2',
    # 'O2*': 'O2',
    # "C1'": 'C1',
    # 'C1*': 'C1',

    "N6": "N6",
    "O2": "O2",
    "C6": "C6",
    "C5": "C5",
    "N7": "N7",
    "C8": "C8",
    "N9": "N9",
    "C4": "C4",
    "N3": "N3",
    "C2": "C2",
    "N1": "N1",
}

ISOGUANINE_ATOM_RES = {
    "N6": 0,
    "O2": 0,
    "C6": 0,
    "C5": 0,
    "N7": 0,
    "C8": 0,
    "N9": 0,
    "C4": 0,
    "N3": 0,
    "C2": 0,
    "N1": 0,
}

ISOGUANINE_REQUIRED_CONDITION = {
    ("N1", "C2", 2.0, 0, 0),
    ("C2", "N3", 2.0, 0, 0),
    ("N3", "C4", 2.0, 0, 0),
    ("C4", "C5", 2.0, 0, 0),
    ("C5", "C6", 2.0, 0, 0),
    ("C6", "N1", 2.0, 0, 0),
    ("C5", "N7", 2.0, 0, 0),
    ("N7", "C8", 2.0, 0, 0),
    ("C8", "N9", 2.0, 0, 0),
    ("N9", "C4", 2.0, 0, 0),
    ("C6", "N6", 2.0, 0, 0),
    ("C2", "O2", 2.0, 0, 0),
}

ISOGUANINE_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': [
        "aC6N1C2",
        "aN1C2N3",
        "aC2N3C4",
        "aN3C4C5",
        "aC4C5C6",
        "aC5C6N1",
        "aN3C4N9",
        "aC6C5N7",
        "aC5C4N9",
        "aC4N9C8",
        "aN9C8N7",
        "aC8N7C5",
        "aN7C5C4",
        "aN6C6N1",
        "aN6C6C5",
        "aO2C2N1",
        "aO2C2N3",
    ]
}

ISOGUANINE_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': []
}

ISOGUANINE_RESTRAINTS = [
    {
        "conditions": [],
        "name": "Base==Isoguanine",
        "restraints": [
            ["dist", "dN1C2", ["N1", "C2"], 1.419, 0.008],
            ["dist", "dC2N3", ["C2", "N3"], 1.351, 0.008],
            ["dist", "dN3C4", ["N3", "C4"], 1.334, 0.008],
            ["dist", "dC4C5", ["C4", "C5"], 1.407, 0.008],
            ["dist", "dC5C6", ["C5", "C6"], 1.401, 0.008],
            ["dist", "dC6N1", ["C6", "N1"], 1.361, 0.008],
            ["dist", "dC5N7", ["C5", "N7"], 1.388, 0.008],
            ["dist", "dN7C8", ["N7", "C8"], 1.312, 0.008],
            ["dist", "dC8N9", ["C8", "N9"], 1.390, 0.008],
            ["dist", "dN9C4", ["N9", "C4"], 1.378, 0.008],
            ["dist", "dC6N6", ["C6", "N6"], 1.337, 0.008],
            ["dist", "dC2O2", ["C2", "O2"], 1.254, 0.008],

            ["angle", "aC6N1C2", ["C6", "N1", "C2"], 126.01, 0.6],
            ["angle", "aN1C2N3", ["N1", "C2", "N3"], 119.57, 0.6],
            ["angle", "aC2N3C4", ["C2", "N3", "C4"], 114.05, 0.6],
            ["angle", "aN3C4C5", ["N3", "C4", "C5"], 129.29, 0.6],
            ["angle", "aC4C5C6", ["C4", "C5", "C6"], 116.54, 0.6],
            ["angle", "aC5C6N1", ["C5", "C6", "N1"], 114.53, 0.6],
            ["angle", "aN3C4N9", ["N3", "C4", "N9"], 126.28, 0.6],
            ["angle", "aC6C5N7", ["C6", "C5", "N7"], 132.04, 0.6],
            ["angle", "aC5C4N9", ["C5", "C4", "N9"], 104.43, 0.6],
            ["angle", "aC4N9C8", ["C4", "N9", "C8"], 106.71, 0.6],
            ["angle", "aN9C8N7", ["N9", "C8", "N7"], 113.39, 0.6],
            ["angle", "aC8N7C5", ["C8", "N7", "C5"], 104.05, 0.6],
            ["angle", "aN7C5C4", ["N7", "C5", "C4"], 111.42, 0.6],
            ["angle", "aN6C6N1", ["N6", "C6", "N1"], 119.78, 0.6],
            ["angle", "aN6C6C5", ["N6", "C6", "C5"], 125.69, 0.6],
            ["angle", "aO2C2N1", ["O2", "C2", "N1"], 116.60, 0.6],
            ["angle", "aO2C2N3", ["O2", "C2", "N3"], 123.83, 0.6],
        ]
    },
]
