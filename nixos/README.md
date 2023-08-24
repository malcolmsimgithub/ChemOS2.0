# NixOS
## Install
The `nixos/chemos` folder contains the `.nix` files needed to install the
software stack presented in our work.

To install them in an existing NixOS system, move the files to your `/etc/nixos`
folder and import the `base.nix` file in your `configuration.nix`:

``` nix
# configuration.nix
imports = [
  ./chemos/base.nix
]
```

And then rebuild your system:

``` sh
nixos-rebuild switch
```

To install Nix or NixOS follow the instructions in the official
[webpage](https://nixos.org/download).
Our configuration has been succesfully tested with NixOS 22.11 and 23.05.
We have decided not to implement our repository as a `flake` as for the moment
of writing this README it still remains as an experimental feature.


## Python

Python packages not present in `nixpkgs` are included as an overlay for
`python3Packages`. The included packages are defined in the `chemos/python.nix`
file. If you need to add additional Python packages to your system, you can
either add them to `scy-python` variable or add the implemented packages in your
own Python definition.

``` nix
# chemos/python.nix
sci-python = pkgs.python3.withPackages (ps: with ps; [
    sila2
    atlas
    olympus
    psycopg
    streamlit
    psycopg2
    # your-package
]);

```

## AiiDA
[AiiDA](https://github.com/aiidateam/aiida-core) has been deployed by using
[mach-nix](https://github.com/DavHau/mach-nix), and uses an independent Python
installation.
To execute an AiiDA workflow use the command `verdi` ([official manual](
https://aiida.readthedocs.io/projects/aiida-core/en/v2.4.0/topics/cli.html)).
If you need to add Python packages to the AiiDA python install include them in
the `requirements` variable found in the `chemos/workflow_manager/aiida.nix` or
execute them by using the desired Python install.

## Database creation
By default, a [PostgreSQL](https://www.postgresql.org/) database named `chemos`
and using `chemos` as both username password will be generated during the
install.
As a security measure, we recommend to change the `ROLE` and `PASSWORD` values
inside `chemos/database/postgresql.nix`:

``` sql
CREATE ROLE <username> WITH LOGIN PASSWORD '<password>' CREATEDB;
CREATE DATABASE chemos;
```

We have decided to not automatically populate the database with our custom
tables, but instead to include the `chemos/database/database.py` Python script
that uses [SQLAlchemy](https://www.sqlalchemy.org/) to prepare the database.
Adapt the script to the needs of your laboratory and then run it with Python
after the install:

``` sh
python chemos/database/database.py
```


# nix-shell
For testing purposes, we added a `nix-shell` containing the packages needed for
the simulations included in this repository. Notice that it does not include the
database neither AiiDA.

To execute it be sure that your computer has a running `nix` installation and
execute it with `nix-shell`:

``` sh
nix-shell chemos/shell.nix
```
