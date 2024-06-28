from PIL import Image
import os

def convert_ppm_to_jpeg(ppm_file_path, jpeg_file_path=None):
    """
    Convert a PPM file to a JPEG file.

    :param ppm_file_path: Path to the input PPM file
    :param jpeg_file_path: Path to save the output JPEG file. If None, save it in the same directory.
    :return: Path to the saved JPEG file
    """
    # Open the PPM file
    with Image.open(ppm_file_path) as img:
        # Convert to RGB mode if necessary
        if img.mode != 'RGB':
            img = img.convert('RGB')

        # Set JPEG file path if not provided
        if jpeg_file_path is None:
            jpeg_file_path = os.path.splitext(ppm_file_path)[0] + '.jpg'

        # Save the image as JPEG
        img.save(jpeg_file_path, 'JPEG')

    return jpeg_file_path

def convert_all_ppm_in_directory(directory_path):
    """
    Convert all PPM files in a given directory to JPEG files.

    :param directory_path: Path to the directory containing PPM files.
    """
    # List all .ppm files in the directory
    ppm_files = [f for f in os.listdir(directory_path) if f.lower().endswith('.ppm')]
    
    # Convert each .ppm file to .jpeg
    for ppm_file in ppm_files:
        ppm_file_path = os.path.join(directory_path, ppm_file)
        jpeg_file_path = os.path.splitext(ppm_file_path)[0] + '.jpg'
        convert_ppm_to_jpeg(ppm_file_path, jpeg_file_path)
        print(f"Converted {ppm_file_path} to {jpeg_file_path}")

# Example usage:
if __name__ == "__main__":
    directory = "./"  # Replace with your directory path
    convert_all_ppm_in_directory(directory)
