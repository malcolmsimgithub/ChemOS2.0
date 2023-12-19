import json
import pickle


def load_json(path_to_file):
    """Loads dictionary (or other json-serializable object) from json file

    Parameters:
        path_to_file (str or pathlib.Path)

    Returns:
        loaded_object (dict or other): loaded object
    """
    with open(path_to_file, 'r') as jsonfile:
        loaded_object = json.load(jsonfile)

    return loaded_object


def save_json(object_to_save, path_to_file):
    """Saves dictionary (or other json-serializable object) to json file

    Parameters:
        object_to_save (dict or other)
        path_to_file (str or pathlib.Path

    Returns:
        None
    """
    with open(path_to_file, 'w') as jsonfile:
        json.dump(object_to_save, jsonfile, indent=2, separators=(",", ": "))


def load_pkl(path_to_file):
    """Loads a python object from a pickle-compressed file.

    Parameters:
        path_to_file (str or pathlib.Path): Path to the target_zone file

    Returns:
        object
    """
    with open(path_to_file, 'rb') as pklfile:
        loaded_object = pickle.load(pklfile)

    return loaded_object


def save_pkl(object_to_save, path_to_file):
    """Saves object to binary pickle file.

    Parameters:
        object_to_save (Any)
        path_to_file (str or pathlib.Path)

    Returns:
        None
    """
    with open(path_to_file, 'wb') as pklfile:
        pickle.dump(object_to_save, pklfile)


def save_csv(object_to_save, path_to_file, header=None):
    """Save dictionary as csv file.

    Parameters:
        object_to_save (dict): Dictionary to save
        path_to_file (Path): Path for target_zone file
        header (list or tuple): column headers to write

    Returns:
        None
    """
    with open(path_to_file, "w") as outfile:
        if header:
            outfile.write(f"{join_any(header, ',')}\n")
        for key in object_to_save:
            if isinstance(object_to_save[key], (str, float, int, bool)):
                outfile.write(f"{key},{object_to_save[key]}\n")
            elif isinstance(object_to_save[key], (list, tuple)):
                outfile.write(f"{key},{join_any(object_to_save[key], ',')}\n")
            elif isinstance(object_to_save[key], dict):
                outfile.write(f"{key},{join_any(object_to_save[key].values(), ',')}\n")


def load_csv(path_to_file, header=False):
    """Load csv file as dictionary (keys: first column; values: all other columns as list).

    Parameters:
        path_to_file (pathlib.Path): Path to the csv file
        header (bool): existance of header in csv file

    Returns:
        file_dict (dict): Dictionary of file contents.
    """
    file_dict = dict
    with open(path_to_file, 'r') as csvfile:
        for i, line in enumerate(csvfile):
            if i == 0:
                if header:
                    continue
            else:
                line_split = line.strip().split(",")
                file_dict[line_split[0]] = [line.split[entry] for entry in range(1,len(line.split))]
    return file_dict


def join_any(iterable, sepsign):
    return f"{sepsign}".join([str(entry) for entry in iterable])
