import logging
from pathlib import Path

def configure_logger(
    logger,
    console_level = logging.WARNING,
    file_level = logging.DEBUG,
    file = None,
) -> logging.Logger:
    """Configures the logging levels and output for a pre-existing logger instance."""
    # Create a standard logging format.
    fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(fmt, datefmt="%Y-%m-%d %H:%M:%S")

    # Setup console logging.
    console_log = logging.StreamHandler()
    console_log.setLevel(console_level)
    console_log.setFormatter(formatter)
    logger.addHandler(console_log)

    # Setup file logging, if a file path is provided
    if file is not None:
        Path(file).parent.mkdir(parents=True, exist_ok=True)
        try:
            file_log = logging.FileHandler(str(file))
        except IsADirectoryError as err:
            raise IsADirectoryError(f"Invalid path '{str(file)}' for file log.") from err
        file_log.setLevel(file_level)
        file_log.setFormatter(formatter)
        logger.addHandler(file_log)
        logger.debug(f"Initialising file log output to {str(file)}.")

    return logger

def setup_logger(
    name = None,
    console_level = logging.WARNING,
    file_level = logging.DEBUG,
    file = None,
):
    """Creates and/or gets a new logger instance 'name' and configures it for use."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    return configure_logger(logger, console_level, file_level, file)
