The model describes pyrolysis intrinsic kinetics. Further physical parameters integration & modelization is expected.
Main program is in main.py

Proposed structure:
/lignin_pyrolysis/:(Root project folder)

├── data/ (contains external data files like reactions.csv and species.csv)
│   ├── reactions.csv
│   └── species.csv
├── modules/ (Contains all the Python modules .py files that make the codebase)
│   ├── __init__.py
│   ├── reactions.py
│   ├── species.py
│   ├── data_handler.py
│   ├── kinetics.py
│   ├── ode_solver.py
│   └── utils.py
├── main.py (The main script that runs the simulation)
└── README.md
