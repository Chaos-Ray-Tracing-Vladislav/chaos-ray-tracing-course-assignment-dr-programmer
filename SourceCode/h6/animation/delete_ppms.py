import os
import glob

def delete_ppm_files():
    # Use glob to find all .ppm files in the current directory
    ppm_files = glob.glob("*.ppm")
    
    # Iterate over the list of .ppm files and delete each one
    for file in ppm_files:
        try:
            os.remove(file)
            print(f"Deleted: {file}")
        except OSError as e:
            print(f"Error deleting {file}: {e}")

if __name__ == "__main__":
    delete_ppm_files()

