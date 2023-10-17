# 🛠️ This file contains general, re-usable Python code that might come in handy from time to time.

import os



def check_if_exist(path):
    '''
    Check if a file at the specified path exists and, if so, prompt the user to overwrite it.

    Parameters:
    - path (str): The path to the file to be checked.

    Returns:
    - bool: True if the file doesn't exist or if the user chooses to overwrite it.
            False if the file exists and the user chooses not to overwrite it.
    '''
    exist = os.path.isfile(path)
    if not exist:
        return True
    else:
        response = input(f'File - {os.path.basename(path)} - already exist. Do you want to overwrite the file (y/n)')
        if response == 'y':
            os.remove(path)
            print(f'File - {os.path.basename(path)} - was removed.')
            return True
        elif response == 'n':
            return False
        else:
            print("Please specify 'y' or 'n'")
            return False
