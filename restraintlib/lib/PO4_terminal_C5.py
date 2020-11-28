
PO4_5_TERMINAL_PDB_CODES = ['A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT', 'DU', 'IC', 'IG']

PO4_5_TERMINAL_ATOM_NAMES = {
    'P': 'P',
    'OP1': 'OP1',
    'O1P': 'OP1',
    'OP2': 'OP2',
    'O2P': 'OP2',
    "OP3": "OP3",
    "O3P": "OP3",
    "O5'": "O5'",
    "O5*": "O5'",
    "C5'": "C5'",
    "C5*": "C5'",
}

PO4_5_TERMINAL_ATOM_RES = {
    'P': 0,
    'OP1': 0,
    'OP2': 0,
    "OP3": 0,
    "O5'": 0,
    "C5'": 0,
}

PO4_5_TERMINAL_REQUIRED_CONDITION = {
    ("P", "O5'", 2.0, 0, 0),
    ("P", "OP1", 2.0, 0, 0),
    ("P", "OP2", 2.0, 0, 0),
    ("P", "OP3", 2.0, 0, 0),
    ("C5'", "O5'", 2.0, 0, 0),
}

PO4_5_TERMINAL_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': ['aO1O2', 'aO1O3', 'aO1O5', 'aO2O3', 'aO2O5', 'aO3O5']
}

PO4_5_TERMINAL_CONDITION_DISTANCE_MEASURE = {
    'measure': 'euclidean_angles',
    'restraint_names': ['tC3O3P4O5', 'tC5O5P4O3']
}

PO4_5_TERMINAL_RESTRAINTS = [
    {
        "conditions": [],
        "name": "PO4_Terminal==C5_0",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 111.7, 1.0],
            ["angle", "aO1O3", ["OP1", "P", "OP3"], 114.0, 0.7],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 107.5, 0.7],
            ["angle", "aO2O3", ["OP2", "P", "OP3"], 112.8, 1.0],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 107.2, 0.8],
            ["angle", "aO3O5", ["OP3", "P", "O5'"], 102.8, 1.2],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.0, 2.2],
            ["dist", "dO1P4", ["OP1", "P"], 1.514, 0.009],
            ["dist", "dO2P4", ["OP2", "P"], 1.52, 0.009],
            ["dist", "dO3P4", ["OP3", "P"], 1.514, 0.01],
            ["dist", "dO5P4", ["O5'", "P"], 1.622, 0.009]
        ]
    },
    {
        "conditions": [],
        "name": "PO4_Terminal==C5_1",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 114.0, 0.7],
            ["angle", "aO1O3", ["OP1", "P", "OP3"], 111.7, 1.0],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 107.5, 0.7],
            ["angle", "aO2O3", ["OP2", "P", "OP3"], 112.8, 1.0],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 102.8, 1.2],
            ["angle", "aO3O5", ["OP3", "P", "O5'"], 107.2, 0.8],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.0, 2.2],
            ["dist", "dO1P4", ["OP1", "P"], 1.514, 0.009],
            ["dist", "dO2P4", ["OP2", "P"], 1.514, 0.01],
            ["dist", "dO3P4", ["OP3", "P"], 1.52, 0.009],
            ["dist", "dO5P4", ["O5'", "P"], 1.622, 0.009]
        ]
    },
    {
        "conditions": [],
        "name": "PO4_Terminal==C5_2",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 111.7, 1.0],
            ["angle", "aO1O3", ["OP1", "P", "OP3"], 112.8, 1.0],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 107.2, 0.8],
            ["angle", "aO2O3", ["OP2", "P", "OP3"], 114.0, 0.7],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 107.5, 0.7],
            ["angle", "aO3O5", ["OP3", "P", "O5'"], 102.8, 1.2],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.0, 2.2],
            ["dist", "dO1P4", ["OP1", "P"], 1.52, 0.009],
            ["dist", "dO2P4", ["OP2", "P"], 1.514, 0.009],
            ["dist", "dO3P4", ["OP3", "P"], 1.514, 0.01],
            ["dist", "dO5P4", ["O5'", "P"], 1.622, 0.009]
        ]
    },
    {
        "conditions": [],
        "name": "PO4_Terminal==C5_3",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 114.0, 0.7],
            ["angle", "aO1O3", ["OP1", "P", "OP3"], 112.8, 1.0],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 102.8, 1.2],
            ["angle", "aO2O3", ["OP2", "P", "OP3"], 111.7, 1.0],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 107.5, 0.7],
            ["angle", "aO3O5", ["OP3", "P", "O5'"], 107.2, 0.8],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.0, 2.2],
            ["dist", "dO1P4", ["OP1", "P"], 1.514, 0.01],
            ["dist", "dO2P4", ["OP2", "P"], 1.514, 0.009],
            ["dist", "dO3P4", ["OP3", "P"], 1.52, 0.009],
            ["dist", "dO5P4", ["O5'", "P"], 1.622, 0.009]
        ]
    },
    {
        "conditions": [],
        "name": "PO4_Terminal==C5_4",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 112.8, 1.0],
            ["angle", "aO1O3", ["OP1", "P", "OP3"], 111.7, 1.0],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 107.2, 0.8],
            ["angle", "aO2O3", ["OP2", "P", "OP3"], 114.0, 0.7],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 102.8, 1.2],
            ["angle", "aO3O5", ["OP3", "P", "O5'"], 107.5, 0.7],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.0, 2.2],
            ["dist", "dO1P4", ["OP1", "P"], 1.52, 0.009],
            ["dist", "dO2P4", ["OP2", "P"], 1.514, 0.01],
            ["dist", "dO3P4", ["OP3", "P"], 1.514, 0.009],
            ["dist", "dO5P4", ["O5'", "P"], 1.622, 0.009]
        ]
    },
    {
        "conditions": [],
        "name": "PO4_Terminal==C5_5",
        "restraints": [
            ["angle", "aO1O2", ["OP1", "P", "OP2"], 112.8, 1.0],
            ["angle", "aO1O3", ["OP1", "P", "OP3"], 114.0, 0.7],
            ["angle", "aO1O5", ["OP1", "P", "O5'"], 102.8, 1.2],
            ["angle", "aO2O3", ["OP2", "P", "OP3"], 111.7, 1.0],
            ["angle", "aO2O5", ["OP2", "P", "O5'"], 107.2, 0.8],
            ["angle", "aO3O5", ["OP3", "P", "O5'"], 107.5, 0.7],
            ["angle", "aP4O5C5", ["P", "O5'", "C5'"], 119.0, 2.2],
            ["dist", "dO1P4", ["OP1", "P"], 1.514, 0.01],
            ["dist", "dO2P4", ["OP2", "P"], 1.52, 0.009],
            ["dist", "dO3P4", ["OP3", "P"], 1.514, 0.009],
            ["dist", "dO5P4", ["O5'", "P"], 1.622, 0.009]
        ]
    }
]

