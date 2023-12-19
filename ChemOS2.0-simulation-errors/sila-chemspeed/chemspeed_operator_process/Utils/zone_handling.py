def merge_zones_to_string(zone_list):
    """Merge individual wells to semicolon-separated zone string.

    Parameters:
        zone_list (list): List of individual wells

    Returns:
        zone_string: Zone string of semicolon-separated wells.
    """
    return ";".join(zone_list)


def split_zones_to_list(zone_string):
    """Split semicolon-separated zone string to individual wells.

    Parameters:
        zone_string (str): Zone string of semicolon-separated wells.

    Returns:
        zone_list (list): List of individual wells
    """
    return zone_string.split(";")


def subtract_zone_strings(zone_string_1, zone_string_2):
    """Perform set-like subtraction of two semicolon-separated zone strings maintaining the order of the original zone string.

    Parameters:
        zone_string_1 (str)
        zone_string_2 (str)

    Returns:
        difference_string (str)
    """

    zone_list_1 = split_zones_to_list(zone_string_1)
    zone_list_2 = split_zones_to_list(zone_string_2)

    difference_list = [well for well in zone_list_1 if well not in zone_list_2]

    return merge_zones_to_string(difference_list)


# TODO: Check if this might be deprecated after re-structuring the whole code