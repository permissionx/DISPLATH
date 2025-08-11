#!/usr/bin/env python3
"""
Material templates for DISPLATH
"""

MATERIALS = {
    "graphene": {
        "name": "Graphene",
        "lattice_constants": {"a": 1.42, "b": 3.35},
        "primary_vectors": "3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b",
        "basis": "0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0",
        "basis_types": "1, 1, 1, 1",
        "elements": {
            1: {"symbol": "C", "ed": 22.0, "displacement_energy": 11.0}
        },
        "pmax": 1.45,
        "vacancy_recover_distance": 1.3,
        "recommended_box_sizes": "10, 20, 10",
        "recommended_lattice_ranges": "0 10; 0 20; 5 6",
        "periodic": "true, true, false",
        "grid_vectors": "a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1"
    },
    "hBN": {
        "name": "Hexagonal Boron Nitride",
        "lattice_constants": {"a": 1.446, "b": 3.35},
        "primary_vectors": "3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b",
        "basis": "0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0",
        "basis_types": "1, 1, 2, 2",
        "elements": {
            1: {"symbol": "B", "ed": 18.3, "displacement_energy": 9.15},
            2: {"symbol": "N", "ed": 18.6, "displacement_energy": 9.3}
        },
        "pmax": 1.45,
        "vacancy_recover_distance": 1.3,
        "recommended_box_sizes": "10, 20, 10",
        "recommended_lattice_ranges": "0 10; 0 20; 5 6",
        "periodic": "true, true, false",
        "grid_vectors": "a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1"
    },
    "MoS2": {
        "name": "Molybdenum Disulfide",
        "lattice_constants": {},
        "primary_vectors": "3.1962223053 0.0 0.0; 0.0 5.5360207558 0.0; 0.0 0.0 23.1298294067",
        "basis": """0.500000000 0.166667000 0.500000000;
         0.000000000 0.666667000 0.500000000;
         0.500000000 0.833333000 0.432137000;
         0.000000000 0.333333000 0.432137000;
         0.500000000 0.833333000 0.567863000;
         0.000000000 0.333333000 0.567863000""",
        "basis_types": "1, 1, 2, 2, 2, 2",
        "elements": {
            1: {"symbol": "S", "ed": 7.8, "displacement_energy": 3.9},
            2: {"symbol": "Mo", "ed": 29.1, "displacement_energy": 14.55}
        },
        "pmax": 1.45,
        "vacancy_recover_distance": 3.0,
        "recommended_box_sizes": "50, 30, 3",
        "recommended_lattice_ranges": "0 50; 0 30; 1 2",
        "periodic": "true, true, false",
        "grid_vectors": "6.0 0.0 0.0; 0.0 6.0 0.0; 0.0 0.0 6.0"
    },
    "Si": {
        "name": "Silicon",
        "lattice_constants": {"a": 5.431},
        "primary_vectors": "a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a",
        "basis": """0.0 0.0 0.0;
         0.5 0.5 0.0;
         0.5 0.0 0.5;
         0.0 0.5 0.5;
         0.25 0.25 0.25;
         0.75 0.75 0.25;
         0.75 0.25 0.75;
         0.25 0.75 0.75""",
        "basis_types": "1, 1, 1, 1, 1, 1, 1, 1",
        "elements": {
            1: {"symbol": "Si", "ed": 20.0, "displacement_energy": 10.0}
        },
        "pmax": "a/2/π^0.5",
        "vacancy_recover_distance": 0.0,
        "recommended_box_sizes": "400, 400, 1005",
        "recommended_lattice_ranges": "0 400; 0 400; 2 1000",
        "periodic": "true, true, true",
        "grid_vectors": "a*1.4 0.0 0.0; 0.0 a*1.4 0.0; 0.0 0.0 a*1.4",
        "temperature": 300.0,
        "debye_temperature": 645.0
    },
    "custom": {
        "name": "Custom Material",
        "lattice_constants": {},
        "primary_vectors": "",
        "basis": "",
        "basis_types": "",
        "elements": {},
        "pmax": 1.45,
        "vacancy_recover_distance": 1.0,
        "recommended_box_sizes": "10, 10, 10",
        "recommended_lattice_ranges": "0 10; 0 10; 0 10",
        "periodic": "true, true, true",
        "grid_vectors": ""
    }
}

# Ion types
IONS = {
    "He": {"symbol": "He", "ed": 0.1, "displacement_energy": 0.1},
    "Ne": {"symbol": "Ne", "ed": 0.1, "displacement_energy": 0.1},
    "Ar": {"symbol": "Ar", "ed": 0.1, "displacement_energy": 0.1},
    "Kr": {"symbol": "Kr", "ed": 0.1, "displacement_energy": 0.1},
    "Xe": {"symbol": "Xe", "ed": 0.1, "displacement_energy": 0.1},
    "B": {"symbol": "B", "ed": 1.0, "displacement_energy": 0.5},
    "C": {"symbol": "C", "ed": 1.0, "displacement_energy": 0.5},
    "N": {"symbol": "N", "ed": 1.0, "displacement_energy": 0.5},
    "O": {"symbol": "O", "ed": 1.0, "displacement_energy": 0.5}
}