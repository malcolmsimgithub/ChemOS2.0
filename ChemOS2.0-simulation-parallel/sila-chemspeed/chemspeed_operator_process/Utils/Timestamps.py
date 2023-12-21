import time


def timestamp_date():
    """Get a timestamp string in the YYYY-MM-DD format.

    Parameters:

    Returns:
        timestamp (str)
    """
    timestamp = time.strftime("%y-%m-%d", time.localtime())
    return timestamp


def timestamp_time():
    """Get a timestamp string in the HH-MM format.

    Parameters:

    Returns:
        timestamp (str)
    """
    timestamp = time.strftime("%H-%M", time.localtime())
    return timestamp


def timestamp_time_precise():
    """Get a timestamp string in the HH-MM-SS format.

    Parameters:

    Returns:
        timestamp (str)
    """
    timestamp = time.strftime("%H-%M-%S", time.localtime())
    return timestamp


def timestamp_datetime():
    """Get a timestamp string in the YYYY-MM-DD_HH-MM format.

    Parameters:

    Returns:
        timestamp (str)
    """
    timestamp = time.strftime(f"{timestamp_date()}_{timestamp_time()}", time.localtime())
    return timestamp
