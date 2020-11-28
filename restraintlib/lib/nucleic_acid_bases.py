
ADENINE_PDB_CODES = ['A', 'DA']

ADENINE_ATOM_NAMES = {
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

    'N9': 'N9',
    'C8': 'C8',
    'N7': 'N7',
    'C5': 'C5',
    'C6': 'C6',
    'N6': 'N6',
    'N1': 'N1',
    'C2': 'C2',
    'N3': 'N3',
    'C4': 'C4',
}

ADENINE_ATOM_RES = {
    'N9': 0,
    'C8': 0,
    'N7': 0,
    'C5': 0,
    'C6': 0,
    'N6': 0,
    'N1': 0,
    'C2': 0,
    'N3': 0,
    'C4': 0,
}

ADENINE_REQUIRED_CONDITION = {
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
}

ADENINE_DISTANCE_MEASURE = {
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
    ]
}

ADENINE_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': []
}

ADENINE_RESTRAINTS = [
    {
        "conditions": [],
        "name": "Base==Adenine",
        "restraints": [
            ["dist", "dN1C2", ["N1", "C2"], 1.339, 0.007],
            ["dist", "dC2N3", ["C2", "N3"], 1.330, 0.007],
            ["dist", "dN3C4", ["N3", "C4"], 1.346, 0.006],
            ["dist", "dC4C5", ["C4", "C5"], 1.382, 0.008],
            ["dist", "dC5C6", ["C5", "C6"], 1.406, 0.008],
            ["dist", "dC6N1", ["C6", "N1"], 1.353, 0.007],
            ["dist", "dC5N7", ["C5", "N7"], 1.388, 0.007],
            ["dist", "dN7C8", ["N7", "C8"], 1.311, 0.007],
            ["dist", "dC8N9", ["C8", "N9"], 1.370, 0.008],
            ["dist", "dN9C4", ["N9", "C4"], 1.374, 0.007],
            ["dist", "dC6N6", ["C6", "N6"], 1.334, 0.007],

            ["angle", "aC6N1C2", ["C6", "N1", "C2"], 118.6, 0.6],
            ["angle", "aN1C2N3", ["N1", "C2", "N3"], 129.4, 0.7],
            ["angle", "aC2N3C4", ["C2", "N3", "C4"], 110.5, 0.6],
            ["angle", "aN3C4C5", ["N3", "C4", "C5"], 126.9, 0.6],
            ["angle", "aC4C5C6", ["C4", "C5", "C6"], 117.1, 0.5],
            ["angle", "aC5C6N1", ["C5", "C6", "N1"], 117.5, 0.5],
            ["angle", "aN3C4N9", ["N3", "C4", "N9"], 127.2, 0.7],
            ["angle", "aC6C5N7", ["C6", "C5", "N7"], 132.2, 0.6],
            ["angle", "aC5C4N9", ["C5", "C4", "N9"], 105.9, 0.4],
            ["angle", "aC4N9C8", ["C4", "N9", "C8"], 105.7, 0.4],
            ["angle", "aN9C8N7", ["N9", "C8", "N7"], 113.9, 0.5],
            ["angle", "aC8N7C5", ["C8", "N7", "C5"], 103.8, 0.4],
            ["angle", "aN7C5C4", ["N7", "C5", "C4"], 110.6, 0.5],
            ["angle", "aN6C6N1", ["N6", "C6", "N1"], 118.6, 0.7],
            ["angle", "aN6C6C5", ["N6", "C6", "C5"], 123.9, 0.7],
        ]
    },
]

GUANINE_PDB_CODES = ['G', 'DG']

GUANINE_ATOM_NAMES = {
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

    "N9": "N9",
    "C8": "C8",
    "N7": "N7",
    "C5": "C5",
    "C6": "C6",
    "O6": "O6",
    "N1": "N1",
    "C2": "C2",
    "N2": "N2",
    "N3": "N3",
    "C4": "C4",
}

GUANINE_ATOM_RES = {
    "N9": 0,
    "C8": 0,
    "N7": 0,
    "C5": 0,
    "C6": 0,
    "O6": 0,
    "N1": 0,
    "C2": 0,
    "N2": 0,
    "N3": 0,
    "C4": 0,
}

GUANINE_REQUIRED_CONDITION = {
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
    ("C6", "O6", 2.0, 0, 0),
    ("C2", "N2", 2.0, 0, 0),
}

GUANINE_DISTANCE_MEASURE = {
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
        "aO6C6N1",
        "aO6C6C5",
        "aN2C2N1",
        "aN2C2N3",
    ]
}

GUANINE_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': []
}

GUANINE_RESTRAINTS = [
    {
        "conditions": [],
        "name": "Base==Guanine",
        "restraints": [
            ["dist", "dN1C2", ["N1", "C2"], 1.372, 0.006],
            ["dist", "dC2N3", ["C2", "N3"], 1.327, 0.005],
            ["dist", "dN3C4", ["N3", "C4"], 1.352, 0.006],
            ["dist", "dC4C5", ["C4", "C5"], 1.379, 0.006],
            ["dist", "dC5C6", ["C5", "C6"], 1.418, 0.008],
            ["dist", "dC6N1", ["C6", "N1"], 1.392, 0.006],
            ["dist", "dC5N7", ["C5", "N7"], 1.388, 0.006],
            ["dist", "dN7C8", ["N7", "C8"], 1.308, 0.006],
            ["dist", "dC8N9", ["C8", "N9"], 1.375, 0.006],
            ["dist", "dN9C4", ["N9", "C4"], 1.374, 0.006],
            ["dist", "dC6O6", ["C6", "O6"], 1.238, 0.007],
            ["dist", "dC2N2", ["C2", "N2"], 1.338, 0.007],

            ["angle", "aC6N1C2", ["C6", "N1", "C2"], 125.4, 0.5],
            ["angle", "aN1C2N3", ["N1", "C2", "N3"], 123.6, 0.5],
            ["angle", "aC2N3C4", ["C2", "N3", "C4"], 112.0, 0.4],
            ["angle", "aN3C4C5", ["N3", "C4", "C5"], 128.6, 0.5],
            ["angle", "aC4C5C6", ["C4", "C5", "C6"], 118.9, 0.4],
            ["angle", "aC5C6N1", ["C5", "C6", "N1"], 111.5, 0.5],
            ["angle", "aN3C4N9", ["N3", "C4", "N9"], 125.8, 0.7],
            ["angle", "aC6C5N7", ["C6", "C5", "N7"], 130.3, 0.5],
            ["angle", "aC5C4N9", ["C5", "C4", "N9"], 105.6, 0.5],
            ["angle", "aC4N9C8", ["C4", "N9", "C8"], 106.2, 0.4],
            ["angle", "aN9C8N7", ["N9", "C8", "N7"], 113.2, 0.4],
            ["angle", "aC8N7C5", ["C8", "N7", "C5"], 104.2, 0.4],
            ["angle", "aN7C5C4", ["N7", "C5", "C4"], 110.8, 0.4],
            ["angle", "aO6C6N1", ["O6", "C6", "N1"], 120.1, 0.5],
            ["angle", "aO6C6C5", ["O6", "C6", "C5"], 128.4, 0.6],
            ["angle", "aN2C2N1", ["N2", "C2", "N1"], 116.5, 0.6],
            ["angle", "aN2C2N3", ["N2", "C2", "N3"], 119.9, 0.6],
        ]
    },
]


URACIL_PDB_CODES = ['U', 'DU']

URACIL_ATOM_NAMES = {
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

    "N1": "N1",
    "C2": "C2",
    "O2": "O2",
    "N3": "N3",
    "C4": "C4",
    "O4": "O4",
    "C5": "C5",
    "C6": "C6",
}

URACIL_ATOM_RES = {
    "N1": 0,
    "C2": 0,
    "O2": 0,
    "N3": 0,
    "C4": 0,
    "O4": 0,
    "C5": 0,
    "C6": 0,
}

URACIL_REQUIRED_CONDITION = {
    ("N1", "C2", 2.0, 0, 0),
    ("C2", "N3", 2.0, 0, 0),
    ("N3", "C4", 2.0, 0, 0),
    ("C4", "C5", 2.0, 0, 0),
    ("C5", "C6", 2.0, 0, 0),
    ("C6", "N1", 2.0, 0, 0),
    ("C2", "O2", 2.0, 0, 0),
    ("C4", "O4", 2.0, 0, 0),
}

URACIL_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': [
        "aC6N1C2",
        "aN1C2N3",
        "aC2N3C4",
        "aN3C4C5",
        "aC4C5C6",
        "aC5C6N1",
        "aO2C2N1",
        "aO2C2N3",
        "aO4C4C5",
        "aO4C4N3",
    ]
}

URACIL_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': []
}

URACIL_RESTRAINTS = [
    {
        "conditions": [],
        "name": "Base==Uracil",
        "restraints": [
            ["dist", "dN1C2", ["N1", "C2"], 1.381, 0.009],
            ["dist", "dC2N3", ["C2", "N3"], 1.373, 0.008],
            ["dist", "dN3C4", ["N3", "C4"], 1.381, 0.008],
            ["dist", "dC4C5", ["C4", "C5"], 1.432, 0.008],
            ["dist", "dC5C6", ["C5", "C6"], 1.337, 0.008],
            ["dist", "dC6N1", ["C6", "N1"], 1.374, 0.008],
            ["dist", "dC2O2", ["C2", "O2"], 1.219, 0.008],
            ["dist", "dC4O4", ["C4", "O4"], 1.231, 0.008],

            ["angle", "aC6N1C2", ["C6", "N1", "C2"], 121.1, 0.5],
            ["angle", "aN1C2N3", ["N1", "C2", "N3"], 114.9, 0.6],
            ["angle", "aC2N3C4", ["C2", "N3", "C4"], 127.0, 0.5],
            ["angle", "aN3C4C5", ["N3", "C4", "C5"], 114.5, 0.6],
            ["angle", "aC4C5C6", ["C4", "C5", "C6"], 119.7, 0.6],
            ["angle", "aC5C6N1", ["C5", "C6", "N1"], 122.7, 0.5],
            ["angle", "aO2C2N1", ["O2", "C2", "N1"], 122.8, 0.7],
            ["angle", "aO2C2N3", ["O2", "C2", "N3"], 122.3, 0.6],
            ["angle", "aO4C4C5", ["O4", "C4", "C5"], 126.0, 0.7],
            ["angle", "aO4C4N3", ["O4", "C4", "N3"], 119.5, 0.7],
        ]
    },
]


THYMINE_PDB_CODES = ['T', 'DT']

THYMINE_ATOM_NAMES = {
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

    "N1": "N1",
    "C2": "C2",
    "O2": "O2",
    "N3": "N3",
    "C4": "C4",
    "O4": "O4",
    "C5": "C5",
    "C5M": "C7",
    "C7": "C7",
    "C6": "C6",
}

THYMINE_ATOM_RES = {
    "N1": 0,
    "C2": 0,
    "O2": 0,
    "N3": 0,
    "C4": 0,
    "O4": 0,
    "C5": 0,
    "C7": 0,
    "C6": 0,
}

THYMINE_REQUIRED_CONDITION = {
    ("N1", "C2", 2.0, 0, 0),
    ("C2", "N3", 2.0, 0, 0),
    ("N3", "C4", 2.0, 0, 0),
    ("C4", "C5", 2.0, 0, 0),
    ("C5", "C6", 2.0, 0, 0),
    ("C6", "N1", 2.0, 0, 0),
    ("C2", "O2", 2.0, 0, 0),
    ("C4", "O4", 2.0, 0, 0),
    ("C7", "C5", 2.0, 0, 0),
}

THYMINE_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': [
        "aC6N1C2",
        "aN1C2N3",
        "aC2N3C4",
        "aN3C4C5",
        "aC4C5C6",
        "aC5C6N1",
        "aO2C2N1",
        "aO2C2N3",
        "aO4C4C5",
        "aO4C4N3",
        "aC7C5C4",
        "aC7C5C6",
    ]
}

THYMINE_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': []
}

THYMINE_RESTRAINTS = [
    {
        "conditions": [],
        "name": "Base==Thymine",
        "restraints": [
            ["dist", "dN1C2", ["N1", "C2"], 1.376, 0.008],
            ["dist", "dC2N3", ["C2", "N3"], 1.372, 0.007],
            ["dist", "dN3C4", ["N3", "C4"], 1.382, 0.008],
            ["dist", "dC4C5", ["C4", "C5"], 1.446, 0.008],
            ["dist", "dC5C6", ["C5", "C6"], 1.340, 0.007],
            ["dist", "dC6N1", ["C6", "N1"], 1.381, 0.007],
            ["dist", "dC2O2", ["C2", "O2"], 1.222, 0.008],
            ["dist", "dC4O4", ["C4", "O4"], 1.229, 0.008],
            ["dist", "dC7C5", ["C7", "C5"], 1.498, 0.006],

            ["angle", "aC6N1C2", ["C6", "N1", "C2"], 121.2, 0.5],
            ["angle", "aN1C2N3", ["N1", "C2", "N3"], 114.7, 0.6],
            ["angle", "aC2N3C4", ["C2", "N3", "C4"], 127.1, 0.5],
            ["angle", "aN3C4C5", ["N3", "C4", "C5"], 115.2, 0.5],
            ["angle", "aC4C5C6", ["C4", "C5", "C6"], 118.1, 0.5],
            ["angle", "aC5C6N1", ["C5", "C6", "N1"], 123.6, 0.5],
            ["angle", "aO2C2N1", ["O2", "C2", "N1"], 123.0, 0.7],
            ["angle", "aO2C2N3", ["O2", "C2", "N3"], 122.3, 0.6],
            ["angle", "aO4C4C5", ["O4", "C4", "C5"], 125.0, 0.7],
            ["angle", "aO4C4N3", ["O4", "C4", "N3"], 119.8, 0.6],
            ["angle", "aC7C5C4", ["C7", "C5", "C4"], 118.7, 0.6],
            ["angle", "aC7C5C6", ["C7", "C5", "C6"], 123.2, 0.6],
        ]
    },
]


CYTOSINE_PDB_CODES = ['C', 'DC']

CYTOSINE_ATOM_NAMES = {
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

    "N1": "N1",
    "C2": "C2",
    "O2": "O2",
    "N3": "N3",
    "C4": "C4",
    "N4": "N4",
    "C5": "C5",
    "C6": "C6",
}

CYTOSINE_ATOM_RES = {
    "N1": 0,
    "C2": 0,
    "O2": 0,
    "N3": 0,
    "C4": 0,
    "N4": 0,
    "C5": 0,
    "C6": 0,
}

CYTOSINE_REQUIRED_CONDITION = {
    ("N1", "C2", 2.0, 0, 0),
    ("C2", "N3", 2.0, 0, 0),
    ("N3", "C4", 2.0, 0, 0),
    ("C4", "C5", 2.0, 0, 0),
    ("C5", "C6", 2.0, 0, 0),
    ("C6", "N1", 2.0, 0, 0),
    ("C2", "O2", 2.0, 0, 0),
    ("C4", "N4", 2.0, 0, 0),
}

CYTOSINE_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': [
        "aC6N1C2",
        "aN1C2N3",
        "aC2N3C4",
        "aN3C4C5",
        "aC4C5C6",
        "aC5C6N1",
        "aO2C2N1",
        "aO2C2N3",
        "aN4C4C5",
        "aN4C4N3",
    ]
}

CYTOSINE_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': []
}

CYTOSINE_RESTRAINTS = [
    {
        "conditions": [],
        "name": "Base==Cytosine",
        "restraints": [
            ["dist", "dN1C2", ["N1", "C2"], 1.395, 0.009],
            ["dist", "dC2N3", ["C2", "N3"], 1.353, 0.007],
            ["dist", "dN3C4", ["N3", "C4"], 1.337, 0.008],
            ["dist", "dC4C5", ["C4", "C5"], 1.424, 0.010],
            ["dist", "dC5C6", ["C5", "C6"], 1.338, 0.008],
            ["dist", "dC6N1", ["C6", "N1"], 1.365, 0.007],
            ["dist", "dC2O2", ["C2", "O2"], 1.240, 0.008],
            ["dist", "dC4N4", ["C4", "N4"], 1.330, 0.008],

            ["angle", "aC6N1C2", ["C6", "N1", "C2"], 120.3, 0.5],
            ["angle", "aN1C2N3", ["N1", "C2", "N3"], 119.1, 0.6],
            ["angle", "aC2N3C4", ["C2", "N3", "C4"], 120.1, 0.5],
            ["angle", "aN3C4C5", ["N3", "C4", "C5"], 121.6, 0.6],
            ["angle", "aC4C5C6", ["C4", "C5", "C6"], 117.5, 0.5],
            ["angle", "aC5C6N1", ["C5", "C6", "N1"], 121.2, 0.6],
            ["angle", "aO2C2N1", ["O2", "C2", "N1"], 118.8, 0.8],
            ["angle", "aO2C2N3", ["O2", "C2", "N3"], 122.0, 0.6],
            ["angle", "aN4C4C5", ["N4", "C4", "C5"], 120.3, 0.7],
            ["angle", "aN4C4N3", ["N4", "C4", "N3"], 118.1, 0.6],
        ]
    },
]
