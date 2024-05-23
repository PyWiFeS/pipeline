import logging
from pathlib import Path

def custom_print(logger, level=logging.INFO):
    def log_print(*args, **kwargs):
        message = " ".join(map(str, args))
        logger.log(level, message)
    return log_print

def configure_logger(
    logger,
    console_level=logging.WARNING,
    file_level=logging.DEBUG,
    file=None,
    format=None,
    datefmt=None,
):
    """Configures the logging levels and output for a pre-existing logger instance."""
    if format is None:
        format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(format, datefmt=datefmt)

    console_log = logging.StreamHandler()
    console_log.setLevel(console_level)
    console_log.setFormatter(formatter)
    logger.addHandler(console_log)

    if file is not None:
        Path(file).parent.mkdir(parents=True, exist_ok=True)
        try:
            file_log = logging.FileHandler(str(file))
        except (IsADirectoryError, PermissionError) as err:
            raise RuntimeError(f"Error creating log file at '{str(file)}'") from err
        file_log.setLevel(file_level)
        file_log.setFormatter(formatter)
        logger.addHandler(file_log)
        logger.debug(f"Initialized file log output to {str(file)}.")

    return logger

def setup_logger(
    name=None,
    console_level=logging.WARNING,
    file_level=logging.INFO,
    file=None,
    format=None,
    datefmt="%Y-%m-%d %H:%M:%S",
):
    """Creates and configures a logger instance for use."""
    if name is None:
        name = 'PyWiFeS'
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    return configure_logger(logger, console_level, file_level, file, format, datefmt)