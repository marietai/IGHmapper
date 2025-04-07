import os

def run_command(command):
    result = os.system(command)
    if result != 0:
        raise Exception(f"Command failed: {command}")
    else:
        print(f"Command executed successfully: {command}")
