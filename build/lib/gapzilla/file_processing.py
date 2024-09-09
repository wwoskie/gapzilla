"""
Functions to manipulate files and directories.
"""

import os


def create_dirs_to_output(path_to_dirs: str | os.PathLike) -> None:
    """
    Create directories to the specified output path if they do not already exist.

    Parameters
    ----------
    path_to_dirs : str or os.PathLike
        The path to the directory or directories that need to be created. This
        can be provided as a string or an `os.PathLike` object.

    Returns
    -------
    None
    """

    if not os.path.isdir(path_to_dirs):
        os.makedirs(path_to_dirs)


def create_output_path(
    input_path: str | os.PathLike,
    output_file_name: str,
    path_to_output_folder: str | os.PathLike,
    suffix_file_name: str = "output",
) -> str | os.PathLike:
    """
    Create an output path based on the input path, output file name, output folder,
    and an optional suffix for the file name.

    Parameters
    ----------
    input_path : str or os.PathLike
        The path of the input file which serves as a basis for deriving the
        output path.
    output_file_name : str
        The name of the output file. If not provided, the name will be derived
        from the `input_path` filename and the `suffix_file_name`.
    path_to_output_folder : str or os.PathLike
        The directory where the output file should be saved. If not provided,
        the directory of the `input_path` will be used. If doesn't exits will be created.
    suffix_file_name : str, optional
        The suffix to append to the file name if `output_file_name` is not
        provided. Default is "output".

    Returns
    -------
    str or os.PathLike
        The complete path where the output file should be saved.
    """

    if path_to_output_folder:
        create_dirs_to_output(path_to_output_folder)
    if not output_file_name and not path_to_output_folder:
        directory, base_name = os.path.split(input_path)
        file_name, file_extension = os.path.splitext(base_name)
        new_file_name = f"{file_name}_{suffix_file_name}{file_extension}"

    elif not output_file_name and path_to_output_folder:
        directory, base_name = os.path.split(input_path)
        directory = os.path.join(directory, path_to_output_folder)
        file_name, file_extension = os.path.splitext(base_name)
        new_file_name = f"{file_name}_{suffix_file_name}{file_extension}"

    elif output_file_name and not path_to_output_folder:
        directory = os.path.split(input_path)[0]
        new_file_name = output_file_name

    elif output_file_name and path_to_output_folder:
        directory = path_to_output_folder
        new_file_name = output_file_name

    output_path = os.path.join(directory, new_file_name)

    return output_path


def modify_first_line(
    input_path: str | os.PathLike, output_path: str | os.PathLike
) -> None:
    """
    Modify the first line of a file if it contains the keyword "LOCUS" and
    save the modified content to a new file. This function fixes PROKKA invalid output.

    Parameters
    ----------
    input_path : str or os.PathLike
        The path to the input file whose content is to be modified.
    output_path : str or os.PathLike
        The path to the output file where the modified content will be saved.

    Returns
    -------
    None

    Notes
    -----
    - This function reads the entire content of the input file into a list of
      lines.
    - If the first line of the file contains the keyword "LOCUS" and the text "bp",
      the function inserts "0 " before "bp" in the first line.
    - The modified content is then written to the output file.
    """

    with open(input_path, "r") as file:
        lines = file.readlines()

    if lines and "LOCUS" in lines[0]:
        first_line = lines[0]
        bp_index = first_line.find("bp")
        if bp_index != -1:
            modified_first_line = first_line[:bp_index] + "0 " + first_line[bp_index:]
            lines[0] = modified_first_line

    with open(output_path, "w") as file:
        file.writelines(lines)


def uniquify_path(path: str | os.PathLike) -> str | os.PathLike:
    """
    Generate a unique file path by appending a counter to the file name if a file with
    the given path already exists.

    Parameters
    ----------
    path : str or os.PathLike
        The initial file path that needs to be uniquified.

    Returns
    -------
    str or os.PathLike
        A unique file path with a counter appended if the original path already exists.
    """

    filename, extension = os.path.splitext(path)
    counter = 1

    while os.path.exists(path):
        path = filename + "_(" + str(counter) + ")" + extension
        counter += 1

    return path
