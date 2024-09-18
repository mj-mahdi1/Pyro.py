The model describes pyrolysis intrinsic kinetics. Further physical parameters integration & modelization is expected.
Main program is in main.py

Proposed structure:
/lignin_pyrolysis/:(Root project folder)

├── data/
│   ├── reactions.csv
│   └── species.csv
├── modules/ 
│   ├── __init__.py
│   ├── reactions.py
│   ├── species.py
│   ├── data_handler.py
│   ├── kinetics.py
│   ├── ode_solver.py
│   └── utils.py
├── main.py 
└── README.md
