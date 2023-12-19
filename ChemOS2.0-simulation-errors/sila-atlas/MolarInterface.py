#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Felix Strieth-Kalthoff
import pandas as pd
import time
import os
from typing import Union, Callable, Tuple, Optional
from pathlib import Path

from molar import Client, ClientConfig
from molar.exceptions import MolarBackendError
import requests
from requests.exceptions import ConnectionError
import magic
from dotenv import load_dotenv


load_dotenv()


class MolarInterface:

    status_possible = ('AVAILABLE', 'ACQUIRED', 'PROCESSING', 'SYNTHESIZED', 'SHIPPED', 'RECEIVED', 'DONE', 'FAILED')

    def __init__(self, db_name: str, fragments: Union[tuple, list, set], db_link_name: str = ""):
        """
        Initializes the client to interact with the MOLAR.

        Args:
            db_name: Name of the database.
            fragments: List or tuple of all database fragments.
            db_link_name: Name of the database in the upload link (not required if no downloading is done)
        """
        self._database = db_name
        self._login()
        self._verify_connection()
        self._fragments: list = list(fragments)
        self._db_link_name = db_link_name

    def _login(self) -> None:
        """
        Creates a user client by logging in to the MOLAR.
        """
        self._config = ClientConfig(
            server_url="https://molar.cs.toronto.edu",
            database_name=self._database,
            email=os.environ["MOLAR_USER"],
            password=os.environ["MOLAR_PASSWORD"]
        )
        self._client = Client(self._config)

    def _verify_connection(self) -> None:
        """
        Verifies connection to the MOLAR by client.test_token().
        Returns if connection was verified. Otherwise rises a value error.
        """
        try:
            test_token = self._client.test_token()
            if test_token["is_active"]:
                return
            else:
                raise KeyError("Problem")

        except (KeyError, MolarBackendError, AttributeError, ConnectionError):
            raise RuntimeError("Database connection could not be properly established!")

    def establish_connection(self, attempt: int = 0) -> None:
        """
        Tries to establish connection to the MOLAR.
        Returns if connection was verified. Otherwise tries again (up to five times with 60 seconds waiting in between).
        Terminates the code otherwise.

        Parameters:
             attempt (int): counter of the attempt (for recursive function calls)
        """
        try:
            self._verify_connection()
            return

        except RuntimeError:
            if attempt > 5:
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print("!!! ERROR – DATABASE COULD NOT BE REACHED AFTER FIVE ATTEMPTS !!!")
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                os._exit(1)
            else:
                print("Problem with the MOLAR connection. Re–trying in one minute. ")
                time.sleep(60)
                self._login()
                attempt = attempt+1
                return self.establish_connection(attempt)

    @staticmethod
    def _run_with_connection(function: Callable) -> Callable:
        """
        Decorator to be applied to the MOLAR query functions.
        Static method, so needs to be applied as @_run_with_connection.__get__(0)

        Verifies MOLAR connection before executing function.

        Catches MolarBackendErrors (e.g. for invalid UUIDs etc.). Returns None in these cases.
        """
        def wrapper(self, *args, **kwargs):
            try:
                self.establish_connection()
                return function(self, *args, **kwargs)
            except MolarBackendError as e:
                print(f"!!! MolarBackendError encountered in function {function.__name__} !!!")
                print(e)
                return None
        return wrapper

    @_run_with_connection.__get__(0)
    def _get_token(self) -> str:
        """
        Get a token for uploading files to the MOLAR.
        """
        return self._client.token

    ################################
    # FETCH DATA FROM THE DATABASE #
    ################################

    @_run_with_connection.__get__(0)
    def get_synthesis_details(self, identifier: str, identifier_type: str = "hid"):
        """
        Fetches the details about a specific synthesis from the MOLAR.

        Parameters:
            identifier (str): Identifier of the target_zone molecule
            identifier_type (str): Which identifier to use ("molecule_id", "hid" (default), "smiles")


        Returns:
            details (pandas.DataFrame): Dataframe of details for the respective experiment.
        """
        request = self._client.query_database(
            types=[
                "synthesis",
                "molecule_molecule.molecule_id",
                "product",
                *self._fragments
            ],
            limit=1000,
            joins=[
                {
                    "type": "molecule_molecule",
                    "on": {
                        "column1": "synthesis.molecule_id",
                        "column2": "molecule_molecule.molecule_id"
                    }
                },
                {
                    "type": "product",
                    "on": {
                        "column1": "product.molecule_id",
                        "column2": "molecule_molecule.molecule_id"
                    }
                },
                *[
                    {
                        "type": frag,
                        "on": {
                            "column1": f"{frag}.molecule_id",
                            "column2": f"molecule_molecule.{frag}"
                        }
                    }
                    for frag in self._fragments
                ]
            ],
            aliases=[
                {
                    "type": "molecule",
                    "alias": "product"
                },
                *[
                    {
                        "type": "molecule",
                        "alias": frag
                    }
                    for frag in self._fragments
                ]
            ],
            filters={
                "type": f"product.{identifier_type}",
                "op": "==",
                "value": identifier
            }
        )

        if request.empty:
            raise KeyError(f"The synthesis {identifier} ({identifier_type}) does not exist in the MOLAR.")
        else:
            return request[
                [
                    "synthesis.synthesis_id",
                    "synthesis.status",
                    "synthesis.molecule_id",
                    *[f"{frag}.{col}" for frag in self._fragments for col in ("hid", "smiles")]
                ]
            ]

    @_run_with_connection.__get__(0)
    def get_syntheses(self, lab: str, instrument: str, status: str) -> Union[None, pd.DataFrame]:
        """
        Fetches the details about all syntheses set to a specific status (e.g. "available") for a specific lab and instrument.

        Parameters:
            lab (str): Name of the lab
            instrument (str): Name of the instrument.
            status (str): Status of the syntheses to be fetched

        Returns:
            pending_syntheses (Union[None, pd.DataFrame])
        """
        request = self._client.query_database(
            types=[
                "synthesis",
                "molecule_molecule.molecule_id",
                "lab.name",
                "machine",
                "product",
                *self._fragments
            ],
            limit=1000,
            joins=[
                {
                    "type": "molecule_molecule",
                    "on": {
                        "column1": "synthesis.molecule_id",
                        "column2": "molecule_molecule.molecule_id"
                    }
                },
                {
                    "type": "lab",
                    "on": {
                        "column1": "synthesis.lab_id",
                        "column2": "lab.lab_id"
                    }
                },
                {
                    "type": "machine",
                    "on": {
                        "column1": "synthesis.machine_id",
                        "column2": "machine.machine_id"
                    }
                },
                {
                    "type": "product",
                    "on": {
                        "column1": "product.molecule_id",
                        "column2": "molecule_molecule.molecule_id"
                    }
                },
                *[
                    {
                        "type": frag,
                        "on": {
                            "column1": f"{frag}.molecule_id",
                            "column2": f"molecule_molecule.{frag}"
                        }
                    }
                    for frag in self._fragments
                ]
            ],
            aliases=[
                {
                    "type": "molecule",
                    "alias": "product"
                },
                *[
                    {
                        "type": "molecule",
                        "alias": frag
                    }
                    for frag in self._fragments
                ]
            ],
            filters={
                "filters": [
                    {
                        "type": "lab.name",
                        "op": "==",
                        "value": lab
                    },
                    {
                        "type": "synthesis.status",
                        "op": "==",
                        "value": status
                    },
                    {
                        "type": "machine.name",
                        "op": "==",
                        "value": instrument
                    }
                ],
                "op": "and"
            }
        )
        if request.empty:
            return None
        else:
            return request[
                [
                    "synthesis.synthesis_id",
                    "synthesis.status",
                    "synthesis.molecule_id",
                    "product.hid",
                    "product.optical_properties",
                    "lab.name",
                    "machine.machine_id",
                    "product.smiles",
                    *[f"{frag}.{col}" for frag in self._fragments for col in ("hid", "smiles")]
                ]
            ]

    @_run_with_connection.__get__(0)
    def get_available_syntheses(self) -> Union[None, pd.DataFrame]:
        """
        Fetches the details about all syntheses set to a specific status (e.g. "available").

        Returns:
            pending_syntheses (Union[None, pd.DataFrame])
        """
        request = self._client.query_database(
            types=[
                "synthesis",
                "molecule_molecule.molecule_id",
                "product",
                *self._fragments
            ],
            limit=1000,
            joins=[
                {
                    "type": "molecule_molecule",
                    "on": {
                        "column1": "synthesis.molecule_id",
                        "column2": "molecule_molecule.molecule_id"
                    }
                },
                {
                    "type": "product",
                    "on": {
                        "column1": "product.molecule_id",
                        "column2": "molecule_molecule.molecule_id"
                    }
                },
                *[
                    {
                        "type": frag,
                        "on": {
                            "column1": f"{frag}.molecule_id",
                            "column2": f"molecule_molecule.{frag}"
                        }
                    }
                    for frag in self._fragments
                ]
            ],
            aliases=[
                {
                    "type": "molecule",
                    "alias": "product"
                },
                *[
                    {
                        "type": "molecule",
                        "alias": frag
                    }
                    for frag in self._fragments
                ]
            ],
            filters={
                "filters": [
                    {
                        "type": "synthesis.status",
                        "op": "==",
                        "value": "AVAILABLE"
                    },
                ],
                "op": "and"
            }
        )

        if request.empty:
            return None
        else:
            return request[
                [
                    "synthesis.synthesis_id",
                    "synthesis.status",
                    "synthesis.molecule_id",
                    "product.hid",
                    "product.smiles",
                    *[f"{frag}.{col}" for frag in self._fragments for col in ("hid", "smiles")]
                ]
            ]

    @_run_with_connection.__get__(0)
    def get_all_syntheses(self) -> Union[None, pd.DataFrame]:
        """
        ONLY FOR DEBUGGING --- THIS FUNCTION DOES NOT HAVE ANY USE IN REAL LIFE!
        """
        request = self._client.query_database(
            types=[
                "synthesis",
                "molecule_molecule.molecule_id",
                "product",
                *self._fragments
            ],
            limit=1000,
            joins=[
                {
                    "type": "molecule_molecule",
                    "on": {
                        "column1": "synthesis.molecule_id",
                        "column2": "molecule_molecule.molecule_id"
                    }
                },
                {
                    "type": "product",
                    "on": {
                        "column1": "product.molecule_id",
                        "column2": "molecule_molecule.molecule_id"
                    }
                },
                *[
                    {
                        "type": frag,
                        "on": {
                            "column1": f"{frag}.molecule_id",
                            "column2": f"molecule_molecule.{frag}"
                        }
                    }
                    for frag in self._fragments
                ]
            ],
            aliases=[
                {
                    "type": "molecule",
                    "alias": "product"
                },
                *[
                    {
                        "type": "molecule",
                        "alias": frag
                    }
                    for frag in self._fragments
                ]
            ]
        )

        if request.empty:
            return None
        else:
            return request[
                [
                    "product.hid",
                    "product.smiles",
                    "product.optical_properties",
                    "synthesis.procedure",
                    *[f"{frag}.{col}" for frag in self._fragments for col in ("hid", "smiles")],
                    "synthesis.synthesis_id",
                    "synthesis.status",
                    "synthesis.molecule_id",
                ]
            ]

    @_run_with_connection.__get__(0)
    def get_molecule(self, identifier: str, identifier_type: str = "hid") -> Union[None, pd.DataFrame]:
        """
        Queries the MOLAR for identifying a specific target_zone molecule.
        Returns the corresponding MOLAR entry.

        Parameters:
            identifier (str): Identifier of the molecule
            identifier_type (str): Which identifier to use ("molecule_id", "hid" (default), "smiles" or "CAS")

        Returns:
            target_mol (Union[None, pd.DataFrame]): Dataframe of the query result
        """
        target_mol = self._client.query_database(
            types="molecule",
            filters={
                "type": f"molecule.{identifier_type}",
                "op": "==",
                "value": identifier
            }
        )

        if target_mol.empty:
            return None
        else:
            return target_mol

    @_run_with_connection.__get__(0)
    def get_all_fragments(self, fragment_type: str) -> Union[None, pd.DataFrame]:
        """
        Queries the MOLAR for identifying all molecules that belong to a certain molecule_type.
        Returns the corresponding MOLAR entries.

        Parameters:
            fragment_type: Name of the molecule_type

        Returns: Dataframe of the query result
        """
        request = self._client.query_database(
            types="molecule",
            joins={
                "type": "molecule_type",
                "on": {
                    "column1": "molecule.molecule_type_id",
                    "column2": "molecule_type.molecule_type_id"
                }
            },
            filters={
                "type": "molecule_type.name",
                "op": "==",
                "value": fragment_type
            },
            limit=1000
        )

        if request.empty:
            return None
        else:
            return 

    @_run_with_connection.__get__(0)
    def get_all_fragments2(self, fragment_type: str) -> Union[None, pd.DataFrame]:
        """
        Queries the MOLAR for identifying all molecules that belong to a certain molecule_type.
        Returns the corresponding MOLAR entries.

        Parameters:
            fragment_type: Name of the molecule_type

        Returns: Dataframe of the query result
        """
        request = self._client.query_database(
            types=[
                "molecule",
                "molecule_type",
                "lab_molecule",
                "lab"
            ],
            limit=1000,
            joins=[
                {
                    "type": "molecule_type",
                    "on": {
                        "column1": "molecule.molecule_type_id",
                        "column2": "molecule_type.molecule_type_id"
                    }
                },
                {
                    "type": "lab_molecule",
                    "on": {
                        "column1": "molecule.molecule_id",
                        "column2": "lab_molecule.molecule_id"
                    }
                },
                {
                    "type": "lab",
                    "on": {
                        "column1": "lab_molecule.lab_id",
                        "column2": "lab.lab_id"
                    }
                }
            ],
            filters={
                "filters": [
                    {
                        "type": "molecule_type.name",
                        "op": "==",
                        "value": fragment_type
                    }
                ],
                "op": "and"
            }
        )

        if request.empty:
            return None
        else:
            return request[[
                "molecule.molecule_id",
                "molecule.hid",
                "molecule.smiles",
                "molecule.commercially_available"
            ]]

    @_run_with_connection.__get__(0)
    def get_available_fragments(self, fragment_type: str, lab: str) -> Union[None, pd.DataFrame]:
        """
        Complex query to get all availabe fragments of a certain fragment type in a certain lab.

        Args:
            fragment_type: Name of the fragment type to be queried (e.g. "fragment_a").
            lab: Name of the lab where the fragment should be available.
        """
        request = self._client.query_database(
            types=[
                "molecule",
                "molecule_type",
                "lab_molecule",
                "lab"
            ],
            limit=1000,
            joins=[
                {
                    "type": "molecule_type",
                    "on": {
                        "column1": "molecule.molecule_type_id",
                        "column2": "molecule_type.molecule_type_id"
                    }
                },
                {
                    "type": "lab_molecule",
                    "on": {
                        "column1": "molecule.molecule_id",
                        "column2": "lab_molecule.molecule_id"
                    }
                },
                {
                    "type": "lab",
                    "on": {
                        "column1": "lab_molecule.lab_id",
                        "column2": "lab.lab_id"
                    }
                }
            ],
            filters={
                "filters": [
                    {
                        "type": "molecule_type.name",
                        "op": "==",
                        "value": fragment_type
                    },
                    {
                        "type": "lab.name",
                        "op": "==",
                        "value": lab
                    },
                    {
                        "type": "lab_molecule.available",
                        "op": "==",
                        "value": True
                    }
                ],
                "op": "and"
            }
        )

        if request.empty:
            return None
        else:
            return request[[
                "molecule.molecule_id",
                "molecule.hid",
                "molecule.smiles",
                "molecule.commercially_available"
            ]]

    @_run_with_connection.__get__(0)
    def _get_uuid(self, table: str, identifier: str, identifier_type: str) -> Union[str, None]:
        """
        Fetches the uuid of a lab/molecule/machine/molecule_type from the MOLAR.

        Parameters:
            table (str): Name of the table to get the uuid from.
            identifier (str): Identifier of the respective entry.
            identifier_type (str): Identifier type (name of the table column).

        Returns:
            uuid (str): UUID of the lab/molecule/machine/molecule_type
        """
        entry = self._client.query_database(
            types=table,
            filters={
                "type": f"{table}.{identifier_type}",
                "op": "==",
                "value": identifier
            }
        )

        if entry.empty:
            return None
        else:
            return entry.at[0, f"{table}_id"]

    @_run_with_connection.__get__(0)
    def _get_synthesis_uuid_from_molecule(self, identifier: str, identifier_type: str = "hid") -> Union[str, None]:
        """
        Fetches the uuid of a synthesis from the corresponding molecule identifier.

        Parameters:
            identifier (str): Identifier of the molecule to get the uuid from.
            identifier_type (str): Which identifier to use ("molecule_id", "hid" (default), "smiles")

        Returns:
            uuid (str): UUID of the synthesis
        """
        synthesis = self._client.query_database(
            types=[
                "synthesis",
                "molecule"
            ],
            joins={
                "type": "molecule",
                "on": {
                    "column1": "molecule.molecule_id",
                    "column2": "synthesis.molecule_id"
                }
            },
            filters={
                "type": f"molecule.{identifier_type}",
                "op": "==",
                "value": identifier
            }
        )

        if synthesis.empty:
            return None
        else:
            return synthesis.at[0, "synthesis.synthesis_id"]

    ##############################
    # WRITE DATA TO THE DATABASE #
    ##############################

    @_run_with_connection.__get__(0)
    def create_target_compound(
            self,
            fragments: Union[tuple, list],
            smiles: str,
            procedure: str,
            identifier_type: str = "hid"
    ) -> Tuple[str, str, str]:
        """
        Creates all entries for the synthesis of a new target_zone compound in the MOLAR.
            - creates the entries in the molecule, molecule_molecule, and synthesis table, linked via the molecule_id
            - entries are linked to the fragments
            - sets the synthesis status to "AVAILABLE"

        Parameters:
            fragments (tuple or list): Collection of all fragment identifiers (e.g. human id, MOLAR id, SMILES)
            smiles (str): SMILES string of the target compound
            procedure (str): Name of the synthesis procedure
            identifier_type (str): Type of identifier of the fragments

        Returns:
            mol (dict): Molar event source_zone of the molecule entry generation
            mol_mol (dict): Molar event source_zone of the molecule_molecule entry generation
            synthesis (dict): Molar event source_zone of the synthesis entry generation
        """
        # Get Fragment Details and Create the Target Molecule Identifier
        fragment_details = {
            fragment: self.get_molecule(fragment, identifier_type=identifier_type)
            for fragment in fragments
        }

        mol_hid: str = "".join(frag_details.at[0, "hid"] for frag_details in fragment_details.values())

        # Check if the Entry for the Target Molecule Exists, else Create it

        mol_uuid = self._get_uuid("molecule", mol_hid, "hid")

        if not mol_uuid:
            mol = self._client.create_entry(
                type="molecule",
                data={
                    "smiles": smiles,
                    "molecule_type_id": self._get_uuid("molecule_type", "abc_molecule", "name"),
                    "optical_properties": {},
                    "commercially_available": False,
                    "hid": "".join(frag_details.at[0, "hid"] for frag_details in fragment_details.values())
                }
            )
            mol_uuid = mol["uuid"]

        # Check if the Entry for the Target Molecule_Molecule Exists, else Create it

        mol_mol_uuid = self._get_uuid("molecule_molecule", mol_uuid, "molecule_id")

        if not mol_mol_uuid:
            mol_mol = self._client.create_entry(
                type="molecule_molecule",
                data={
                    "molecule_id": mol_uuid,
                    **{
                        frag: frag_details.at[0, "molecule_id"]
                        for frag, frag_details in zip(self._fragments, fragment_details.values())
                    }
                }
            )
            mol_mol_uuid = mol_mol["uuid"]

        # Create the Entry for the Target Synthesis

        synthesis = self._client.create_entry(
            type="synthesis",
            data={
                "molecule_id": mol_uuid,
                "status": "AVAILABLE",
                "procedure": procedure
            }
        )

        return mol_uuid, mol_mol_uuid, synthesis["uuid"]

    @_run_with_connection.__get__(0)
    def update_synthesis_status(self, identifier: str, status: str, identifier_type: str = "hid") -> None:
        """
        Updates the status of a synthesis run.

        Parameters:
            identifier (str): Identifier of the molecule
            status (str): New status to be uploaded. P
            identifier_type (str): Which identifier to use ("molecule_id", "hid" (default), "smiles")
        """
        assert status in self.status_possible

        synthesis_id = self._get_synthesis_uuid_from_molecule(identifier, identifier_type)

        self._client.update_entry(
            uuid=synthesis_id,
            type="synthesis",
            data={"status": status}
        )

    @_run_with_connection.__get__(0)
    def claim_synthesis(self, identifier: str, instrument: str, lab: str, status: str = "ACQUIRED", identifier_type: str = "hid"):
        """
        Assigns a synthesis run to a specific lab and instrument.

        Parameters:
            identifier (str): Identifier of the molecule
            instrument (str): Name of the instrument to be included (currently: "The Machine", "ChemSpeed", "Chemputer")
            lab (str): Name of the lab to be included (currently: "Toronto", "Illinois", "Vancouver", "Glasgow", if respective instrument is available)
            status (str): New status to be uploaded.
            identifier_type (str): Which identifier to use ("molecule_id", "hid" (default), "smiles")
        """
        assert status in self.status_possible

        synthesis_id = self._get_synthesis_uuid_from_molecule(identifier, identifier_type)
        lab_id = self._get_uuid("lab", lab, "name")
        instrument_id = self._get_uuid("machine", instrument, "name")

        self._client.update_entry(
            uuid=synthesis_id,
            type="synthesis",
            data={"status": status,
                  "lab_id": lab_id,
                  "machine_id": instrument_id}
        )

    @_run_with_connection.__get__(0)
    def update_optics_data(self, identifier: str, optics_data: dict, identifier_type: str = "hid") -> None:
        """Uploads the optical characterization data after completion of a synthesis run.
            - Checks for the uuid of the target_zone molecule in the molecule table (molecule.molecule_id).
            - If the target_zone molecule is not available, it attempts to create this entry from the HID.
            - Writes the characterization data to the entry in the molecule table (molecule.optical_properties).
            - Updates the Synthesis status to "DONE" or "FAILED"

        Parameters:
            identifier (str): Identifier of the synthesis run.
            optics_data (dict): Optical properties to be uploaded
            identifier_type (str): Which identifier to use ("molecule_id", "hid" (default), "smiles")
        """
        target_mol = self.get_molecule(identifier, identifier_type)

        if target_mol is not None:
            target_uuid = target_mol.at[0, "molecule_id"]

        elif identifier_type == "hid":  # ATTENTION: CURRENTLY SPECIAL FOR A0001B001C001 HID TYPE
            try:
                new_entry = self.create_target_compound(
                    fragments=(identifier[0:4], identifier[4:8], identifier[8:12]),
                    instrument_id="f158c26e-a94e-4994-9a95-de324aa1da23",  # ATTENTION: ASSIGNS ALL NEW SYNTHESES TO CHEMSPEED IN TORONTO
                    lab_id="08cc75ca-8d1e-4ead-9b56-ad1b6ba83a41"
                )

                target_uuid = new_entry[0]["uuid"]

            except (ValueError, KeyError, MolarBackendError):
                raise KeyError(f"Compound {identifier} does not exist in MOLAR and entry could not be created.")

        else:
            raise KeyError(f"Compound {identifier} does not exist in MOLAR and entry could not be created from this identifier.")


        if optics_data["validation_status"] is True:
            target_status = "DONE"
        else:
            target_status = "FAILED"

        self._client.update_entry(
            uuid=target_uuid,
            type="molecule",
            data={
                "optical_properties": optics_data,
            }
        )

        self.update_synthesis_status(identifier, target_status, identifier_type)

    @_run_with_connection.__get__(0)
    def _send_file(self, file_path: Path, synthesis_uuid: str) -> None:
        """
        --- Method from Theo (only slight modifications) ---
        Sends a file from a file location (given by file_path) to the MOLAR.
        Links it to the synthesis given by synthesis_uuid.
        Raises an Exception if anything went wrong during file upload (status code not 200).

        Parameters:
            file_path (Path): Path to the file to be uploaded
            synthesis_uuid (str): UUID of the synthesis to be linked to
        """
        token = self._get_token()

        out = requests.post(
            f"https://molar.cs.toronto.edu/{self._db_link_name}/v1/upload-experiment",
            params={"synthesis_id": synthesis_uuid},
            headers={"Authorization": f"Bearer {token}"},
            files={
                "file": (
                    file_path.name,
                    open(file_path, "rb"),
                    magic.from_file(str(file_path), mime=True),
                )
            },
        )
        if out.status_code == 200:
            print(f"Successfully uploaded {file_path}.")
            return
        raise MolarBackendError(out.status_code, f"Something went wrong: {out.text}")

    @_run_with_connection.__get__(0)
    def upload_hplc_data(self, data_archive: Path, identifier: str, identifier_type: str = "hid"):
        """
        Uploads an HPLC data archive to the MOLAR and links it to a synthesis:
            - gets the synthesis uuid
            - uploads the file to the MOLAR using self._send_file()

        Parameters:
            data_archive (Path): Path to the HPLC-MS data archive
            identifier (str): Identifier of the target_zone compound to be synthesized
            identifier_type (str): Which identifier to use ("molecule_id", "hid" (default), "smiles")
        """
        synthesis_uuid = self._get_synthesis_uuid_from_molecule(identifier, identifier_type)
        if synthesis_uuid is None:
            raise ValueError(f"The molecule entry corresponding to {identifier} in the MOLAR does not exist.")
        self._send_file(data_archive, synthesis_uuid)

    ####################################
    # DELETE ENTRIES FROM THE DATABASE #
    ####################################

    def delete_entry(self, table: str, identifier: str, identifier_type: str = "hid") -> None:
        """
        Deletes an entry from the database.

        Args:
            table: Name of the table.
            identifier: Identifier of the entry to be deleted.
            identifier_type: Type of the identifier (column name in the target table).
        """
        uuid = self._get_uuid(table, identifier, identifier_type)

        self._client.delete_entry(
            type=table,
            uuid=uuid
        )

    def delete_target_compound(self, identifier: str, identifier_type: str = "hid") -> None:
        """
        Deletes all entries related to a target compound (entries in the molecule, molecule_molecule, and synthesis
        tables) from the database.

        Args:
            identifier: Identifier of the target molecule.
            identifier_type: Type of the identifier (column name in the target table).
        """
        mol_uuid = self._get_uuid("molecule", identifier, identifier_type)
        mol_mol_uuid = self._get_uuid("molecule_molecule", mol_uuid, "molecule_id")
        synthesis_uuid = self._get_uuid("synthesis", mol_uuid, "molecule_id")

        self.delete_entry("synthesis", synthesis_uuid, "synthesis_id")
        self.delete_entry("molecule_molecule", mol_mol_uuid, "molecule_molecule_id")
        self.delete_entry("molecule", mol_uuid, "molecule_id")