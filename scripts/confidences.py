import glob
import pandas as pd
import numpy as np
import argparse


def get_confidences(pdb_files: list) -> list:
    """
    Read AlphaFold PDB file contents to find confidence info
    @param: pdbs glob list of AlphaFold PDB files
    @return: confidence_info list containing the dictionary of confidence info
    """
    counter = 1
    confidence_info = []
    for pdb_file in pdb_files:
        print(f"At {counter}/{len(pdb_files)}")
        file = open(pdb_file)
        lines = file.readlines()
        file.close()
        confidence_scores = []
        confidence_dict = {}
        for line in lines:
            if line.startswith("ATOM"):
                current_line = line.split()
                confidence_score = float(current_line[-2])
                confidence_scores.append(confidence_score)
        confidence_dict["mean"] = np.mean(confidence_scores)
        confidence_dict["std"] = np.std(confidence_scores)
        confidence_dict["filename"] = pdb_file
        confidence_info.append(confidence_dict)
        counter += 1
    return confidence_info
                                          
                                      
def filter_confidences(pdb_files: list, low_confidence: bool) -> pd.DataFrame:
    """
    Filter out AlphaFold structures with less than 90% confidence
    @param pdb_files: glob list of AlphaFold PDBs
    @param low_confidence: True: get confidences below 90.0, False: get confidences above 90.0
    @return good_confidence_df: dataframe containing PDBs with either low or high confidence
    """
    confidence_list = get_confidences(pdb_files)
    confidence_dataframe = pd.DataFrame.from_dict(confidence_list)
    if low_confidence: 
        confidence_filter = confidence_dataframe["mean"] <= 90.0
    elif not low_confidence:
        confidence_filter = confidence_dataframe["mean"] > 90.0
    return confidence_dataframe.where(confidence_filter).dropna()


def check_positive(input) -> float:
    """
    Type-checking for per-residue confidence scores
    @param input: user input from command-line arguments
    @return: checked per-residue confidence score
    """
    try:
        confidence = float(input)
        if confidence <= 0:
            raise argparse.ArgumentTypeError(f"{confidence} is an invalid per-residue confidence score")
    except ValueError:
        print("Error: per-residue confidence score should be a number")
    except argparse.ArgumentTypeError as message:
        print(message)
    return confidence


def commandline_arguments() -> argparse.Namespace:
    """
    Parser for commandline arguments
    @return: commandline arguments from user
    """
    parser = argparse.ArgumentParser(description="filter AlphaFold 2 structures based on per-residue confidence scores")
    parser.add_argument("-l",
                        "--low",
                        help="get low-confidence AlphaFold 2 structures, default behaviour is to only get structures with confidence above 90",
                        action="store_true")
    return parser.parse_args()


def main():
    arguments = commandline_arguments()
    low_confidence = arguments.low
    pdbs = glob.glob("../data/alphafold/pdb_files/*.pdb")
    confidence_df = filter_confidences(pdbs, low_confidence)
    if low_confidence:
        savename = "low_confidences.csv"
    else:
        savename = "confidences.csv"
    confidence_df.to_csv(f"../data/alphafold/{savename}")


if __name__ == "__main__":
    main()
