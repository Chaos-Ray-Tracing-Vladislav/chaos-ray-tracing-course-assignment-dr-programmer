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

# Example usage:
if __name__ == "__main__":
    for i in range(120):
        ppm_file = "task5_animation_frame_" + str(i+1) + "_1920_1080.ppm"  # Replace with your PPM file path
        jpeg_file = "task5_animation_frame_" + str(i+1) + "_1920_1080.jpeg"  # Optionally replace with desired JPEG file path
        output_path = convert_ppm_to_jpeg(ppm_file, jpeg_file)
        print(f"Converted PPM to JPEG: {output_path}")

