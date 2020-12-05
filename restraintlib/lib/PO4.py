
PO4_PDB_CODES = ['A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT', 'DU', 'IC', 'IG']

PO4_ATOM_NAMES = {
    'P': 'P',
    'OP1': 'OP1',
    'O1P': 'OP1',
    'OP2': 'OP2',
    'O2P': 'OP2',
    "O5'": "O5'",
    "O5*": "O5'",
    "O3'": "O3'",
    "O3*": "O3'",
    "C5'": "C5'",
    "C5*": "C5'",
    "C3'": "C3'",
    "C3*": "C3'",
}

PO4_ATOM_RES = {
    'P': 0,
    'OP1': 0,
    'OP2': 0,
    "O5'": 0,
    "O3'": -1,
    "C5'": 0,
    "C3'": -1,
}

PO4_REQUIRED_CONDITION = {
    ("P", "O3'", 2.0, 0, -1),
    ("P", "O5'", 2.0, 0, 0),
    ("P", "OP1", 2.0, 0, 0),
    ("P", "OP2", 2.0, 0, 0),
    ("C3'", "O3'", 2.0, -1, -1),
    ("C5'", "O5'", 2.0, 0, 0),
}

PO4_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': ['aO1O2', 'aO1O3', 'aO1O5', 'aO2O3', 'aO2O5', 'aO3O5']
}

PO4_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': ['tC3O3P4O5', 'tC5O5P4O3']
}

PO4_RESTRAINTS = [
    {
        "conditions": [
            ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], -66.636, 7.779],
            ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], 171.37, 14.971]
        ],
        "name": "PO4==AA_0",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 117.6, 1.2],
            ["angle", "aO1O3", ["OP1", "P", "O3'"], 106.2, 1.1],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 110.2, 1.3],
            ["angle", "aO2O3", ["OP2", "P", "O3'"], 112.2, 1.0],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 109.3, 0.9],
            ["angle", "aO3O5", ["O3'", "P", "O5'"], 99.9, 0.7],
            ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.2, 1.5],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 121.7, 3.0],
            ["dist", "dO1P4", ["OP1", "P"], 1.487, 0.01],
            ["dist", "dO2P4", ["OP2", "P"], 1.483, 0.01],
            ["dist", "dO3P4", ["O3'", "P"], 1.601, 0.008],
            ["dist", "dO5P4", ["O5'", "P"], 1.591, 0.004],
            ["dist", "dO3C3", ["O3'", "C3'"], 1.422, 0.010],
            ["dist", "dO5C5", ["O5'", "C5'"], 1.428, 0.013]
        ]
    },
    {
        "conditions": [
            ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], 171.37, 14.971],
            ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], -66.636, 7.779]
        ],
        "name": "PO4==AA_1",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 117.6, 1.2],
            ["angle", "aO1O3", ["OP1", "P", "O3'"], 109.3, 0.9],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 112.2, 1.0],
            ["angle", "aO2O3", ["OP2", "P", "O3'"], 110.2, 1.3],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 106.2, 1.1],
            ["angle", "aO3O5", ["O3'", "P", "O5'"], 99.9, 0.7],
            ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.2, 1.5],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 121.7, 3.0],
            ["dist", "dO1P4", ["OP1", "P"], 1.483, 0.01],
            ["dist", "dO2P4", ["OP2", "P"], 1.487, 0.01],
            ["dist", "dO3P4", ["O3'", "P"], 1.601, 0.008],
            ["dist", "dO5P4", ["O5'", "P"], 1.591, 0.004],
            ["dist", "dO3C3", ["O3'", "C3'"], 1.422, 0.010],
            ["dist", "dO5C5", ["O5'", "C5'"], 1.428, 0.013]
        ]
    },
    {
        "conditions": [
            ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], -171.37, 14.971],
            ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], 66.636, 7.779]
        ],
        "name": "PO4==AA_2",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 117.6, 1.2],
            ["angle", "aO1O3", ["OP1", "P", "O3'"], 110.2, 1.3],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 106.2, 1.1],
            ["angle", "aO2O3", ["OP2", "P", "O3'"], 109.3, 0.9],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 112.2, 1.0],
            ["angle", "aO3O5", ["O3'", "P", "O5'"], 99.9, 0.7],
            ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.2, 1.5],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 121.7, 3.0],
            ["dist", "dO1P4", ["OP1", "P"], 1.487, 0.01],
            ["dist", "dO2P4", ["OP2", "P"], 1.483, 0.01],
            ["dist", "dO3P4", ["O3'", "P"], 1.601, 0.008],
            ["dist", "dO5P4", ["O5'", "P"], 1.591, 0.004],
            ["dist", "dO3C3", ["O3'", "C3'"], 1.422, 0.010],
            ["dist", "dO5C5", ["O5'", "C5'"], 1.428, 0.013]
        ]
    },
    {
        "conditions": [
            ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], 66.636, 7.779],
            ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], -171.37, 14.971]
        ],
        "name": "PO4==AA_3",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 117.6, 1.2],
            ["angle", "aO1O3", ["OP1", "P", "O3'"], 112.2, 1.0],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 109.3, 0.9],
            ["angle", "aO2O3", ["OP2", "P", "O3'"], 106.2, 1.1],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 110.2, 1.3],
            ["angle", "aO3O5", ["O3'", "P", "O5'"], 99.9, 0.7],
            ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.2, 1.5],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 121.7, 3.0],
            ["dist", "dO1P4", ["OP1", "P"], 1.483, 0.01],
            ["dist", "dO2P4", ["OP2", "P"], 1.487, 0.01],
            ["dist", "dO3P4", ["O3'", "P"], 1.601, 0.008],
            ["dist", "dO5P4", ["O5'", "P"], 1.591, 0.004],
            ["dist", "dO3C3", ["O3'", "C3'"], 1.422, 0.010],
            ["dist", "dO5C5", ["O5'", "C5'"], 1.428, 0.013]
        ]
    },
    {
        "conditions": [
            ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], -69.896, 9.625],
            ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], -68.72, 8.686]
        ],
        "name": "PO4==AS_0",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 119.9, 1.6],
            ["angle", "aO1O3", ["OP1", "P", "O3'"], 104.5, 0.9],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 110.3, 0.8],
            ["angle", "aO2O3", ["OP2", "P", "O3'"], 111.5, 1.1],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 105.2, 0.8],
            ["angle", "aO3O5", ["O3'", "P", "O5'"], 104.2, 1.5],
            ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.7, 2.9],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.3, 1.5],
            ["dist", "dO1P4", ["OP1", "P"], 1.484, 0.012],
            ["dist", "dO2P4", ["OP2", "P"], 1.478, 0.01],
            ["dist", "dO3P4", ["O3'", "P"], 1.603, 0.014],
            ["dist", "dO5P4", ["O5'", "P"], 1.594, 0.009],
            ["dist", "dO3C3", ["O3'", "C3'"], 1.438, 0.007],
            ["dist", "dO5C5", ["O5'", "C5'"], 1.437, 0.017]
        ]
    },
    {
        "conditions": [
            ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], -68.72, 8.686],
            ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], -69.896, 9.625]
        ],
        "name": "PO4==AS_1",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 119.9, 1.6],
            ["angle", "aO1O3", ["OP1", "P", "O3'"], 105.2, 0.8],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 111.5, 1.1],
            ["angle", "aO2O3", ["OP2", "P", "O3'"], 110.3, 0.8],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 104.5, 0.9],
            ["angle", "aO3O5", ["O3'", "P", "O5'"], 104.2, 1.5],
            ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.7, 2.9],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.3, 1.5],
            ["dist", "dO1P4", ["OP1", "P"], 1.478, 0.01],
            ["dist", "dO2P4", ["OP2", "P"], 1.484, 0.012],
            ["dist", "dO3P4", ["O3'", "P"], 1.603, 0.014],
            ["dist", "dO5P4", ["O5'", "P"], 1.594, 0.009],
            ["dist", "dO3C3", ["O3'", "C3'"], 1.438, 0.007],
            ["dist", "dO5C5", ["O5'", "C5'"], 1.437, 0.017]
        ]
    },
    {
        "conditions": [
            ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], 68.72, 8.686],
            ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], 69.896, 9.625]
        ],
        "name": "PO4==AS_2",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 119.9, 1.6],
            ["angle", "aO1O3", ["OP1", "P", "O3'"], 110.3, 0.8],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 104.5, 0.9],
            ["angle", "aO2O3", ["OP2", "P", "O3'"], 105.2, 0.8],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 111.5, 1.1],
            ["angle", "aO3O5", ["O3'", "P", "O5'"], 104.2, 1.5],
            ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.7, 2.9],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.3, 1.5],
            ["dist", "dO1P4", ["OP1", "P"], 1.484, 0.012],
            ["dist", "dO2P4", ["OP2", "P"], 1.478, 0.01],
            ["dist", "dO3P4", ["O3'", "P"], 1.603, 0.014],
            ["dist", "dO5P4", ["O5'", "P"], 1.594, 0.009],
            ["dist", "dO3C3", ["O3'", "C3'"], 1.438, 0.007],
            ["dist", "dO5C5", ["O5'", "C5'"], 1.437, 0.017]
        ]
    },
    {
        "conditions": [
            ["torsion", "tC3O3P4O5", ["O5'", "P", "O3'", "C3'"], 69.896, 9.625],
            ["torsion", "tC5O5P4O3", ["O3'", "P", "O5'", "C5'"], 68.72, 8.686]
        ],
        "name": "PO4==AS_3",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 119.9, 1.6],
            ["angle", "aO1O3", ["OP1", "P", "O3'"], 111.5, 1.1],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 105.2, 0.8],
            ["angle", "aO2O3", ["OP2", "P", "O3'"], 104.5, 0.9],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 110.3, 0.8],
            ["angle", "aO3O5", ["O3'", "P", "O5'"], 104.2, 1.5],
            ["angle", "aP4O3C3", ["P", "O3'", "C3'"], 120.7, 2.9],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.3, 1.5],
            ["dist", "dO1P4", ["OP1", "P"], 1.478, 0.01],
            ["dist", "dO2P4", ["OP2", "P"], 1.484, 0.012],
            ["dist", "dO3P4", ["O3'", "P"], 1.603, 0.014],
            ["dist", "dO5P4", ["O5'", "P"], 1.594, 0.009],
            ["dist", "dO3C3", ["O3'", "C3'"], 1.438, 0.007],
            ["dist", "dO5C5", ["O5'", "C5'"], 1.437, 0.017]
        ]
    }
]