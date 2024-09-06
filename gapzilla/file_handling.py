import os


def create_dirs_to_output(path_to_dirs):
    if not os.path.isdir(path_to_dirs):
        os.makedirs(path_to_dirs)


def create_output_path(
    input_path, output_file_name, path_to_output_folder, suffix_file_name="output"
):
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


def modify_first_line(input_path, output_path):
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
