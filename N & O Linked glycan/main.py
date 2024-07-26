import argparse
import backend
import logging
import pandas as pd

def setup_logger(log_path: str, name: str = "logger") -> logging.Logger:
    """Configures the root logger.
    Parameters
    ----------
    log_path : str
        The filepath to the log handler.
    name : str (default: "logger")
        The name of the root logger.
    Returns
    -------
    logging.Logger
        The root logger.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename=log_path, encoding="utf-8", mode="w")
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def main():
    parser = argparse.ArgumentParser(prog = "command line tool")
    parser.add_argument("-n", "--nocache", action= "store_false", help= "use the cache files or re-pull from the API")

    options= parser.parse_args()

    logger = setup_logger("./logfile.log")

    backend.all_supersearch_list_data(options.nocache, logger)
    backend.get_all_protein_data(options.nocache, logger)
    backend.get_glycan_data(logger)
    df = pd.read_csv("./output/supersearch_results.tsv", delimiter="\t")
    backend.sites_data(logger, df)



if __name__ == "__main__":
    main()
