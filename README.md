# ChemOS 2.0

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


`ChemOS 2.0` is an orchestration framework for automating chemistry labs. 

This repository contains code used in the creation of ChemOS 2.0, as well as data collected during development.
https://chemrxiv.org/engage/chemrxiv/article-details/64cbe80adfabaf06ffa61204

There is also a demontstrative collection of simulators to showcase the functioning of ChemOS 2.0.

Output data for the machine's runs can be found in the folder "the-machine-data".

Output data for the closed loop can be found in the folder "closed-loop-data".

ChemOS 2.0 is currently under peer review. Additional changes to this repository's code following publication

NOTE: Specific implementations of ChemOS 2.0 require a unique customization of the architecture for a specific experimental setting. 
      repository does not contain a complete package for immediate deployment


# Contact

For all enquiries about the code, please contact [Malcolm] as malcolm.sim@mail.utoronto.ca
 of the sila2 servers used by chemos to conduct a closed-loop experimental campaign.

This code is for demonstration purposes only and is meant to showcase the functionality of ChemOS 2.0

## Installation of simulator code

Install The simulation demo code from source by either cloning this repository or downloading it as a zip file.
All code for the simulations may be found in the folder ChemOS2.0-simulation

Once the repository is downloaded, you must use the nixos configuration, or open a nix shell as instructed in the README for https://github.com/malcolmsimgithub/ChemOS2.0/tree/master/nixos. You can also find instructions for installation of the database here.


IMPORTANT: this code assumes that postgresql is running on port 5432 ((default))

Once Postgresql is installed and you nix environment is configured, please make sure to fill out the script ChemOS2.0-simulation/dblogin.py as well as ChemOS2.0-simulation/streamlit/dblogin.py with your database login credentials.
The default user and password for this code is "chemos"

```Python

def get_login():
    dbname = "chemos"
    dbuser = "chemos"
    dbpassword = "chemos"
    return dbname, dbuser, dbpassword

```


## Usage of simulators


IMPORTANT:  the simulators are primarily meant to showcase the functionality of SiLA2 as an interface for laboratory hardware. They are not meant for deployment or any
real-world usage. The simulators allow users to send commands to virtual clones of the hardware used in the organic-laser campaign of the ChemOS 2.0 project. like their deployment counterparts, they utilize SiLA2 to receive and execute commands from ChemOS2.0, and stream output data back to ChemOS2.0. Whilst the SiLA2 servers for the simulators are all hosted locally on the same device, note that for the actual deployment, each SiLA2 server was hosted on a separate computer connected to its respective hardware. The simulators return (identical) sample data from each instrument, which is automatically logged in the database (which should be configured according to the instructions above).

Next, we can start of the sila2 servers:

```Bash
cd ChemOS2.0-simulation
python start_sila_servers.py
```

To start the gui:

```Bash
cd streamlit
streamlit run Hello.py
```


Finally, to run the "closed loop" using all of the simulators:

```Bash
python laserworkflow.py
```

One can look at results in the database using psql. Data generated is dummy data and is identical for all parameters.

example job files to use with the GUI can be found in the folder ChemOS2.0-simulation/job_files

if a sila2 server crashes, one can restart the servers using the command once more. this will usually clear all of the folders and solve runtime errors.
```Bash
python start_sila_servers.py
```

![Screenshot 2023-10-30 at 11 49 53 AM](https://github.com/malcolmsimgithub/ChemOS2.0/assets/110689046/b7121ecf-e29b-4f2e-a91f-ed461dd84be6)

Web application for ChemOS 2.0. A:sidebar for toggling between instruments/features. B:job submission for the HPLC-MS. C:data visualiza- tion of the HLPC-MS’s job results.

# Faulty simulators
The folder ChemOS2.0-simulation-errors contains code for modified HPLC and optics table simulators that can crash randomly. 
It aslo contains some scripts that demonstrate error handling in workflows
These are intended to demonstrate the error handling capabilities of ChemOS 2.0

# Parallel simulators
The folder ChemOS2.0-simulation-parallel contains code for two HPLC simulators, and a script designed to showcase the job parallelization capabilities of ChemOS 2.0 using python's asyncio library.

## echem data
Note that we have included a folder called "echem_data". This code is a part of a thesis that incorporates this project, and is not directly related to 
the ChemOS 2.0 white paper.

## License

Distributed under the [MIT](https://choosealicense.com/licenses/mit/)
 license. See `LICENSE` for more information.

## Contact

Please reach out to [Malcolm](malcolm.sim@mail.utoronto.ca) by email if you have questions.



