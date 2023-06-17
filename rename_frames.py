import os 
import argparse

parser = argparse.ArgumentParser()
# Make --digit_count a positional argument
parser.add_argument('digit_count', type=int, default=3)

args = parser.parse_args()

FRAMES_DIR = './frames'
digit_count = args.digit_count
for frame_name in os.listdir(FRAMES_DIR):
    # Structure of frame_name: outx.png
    # extract x:
    frame_number = frame_name[3:-4]
    frame_number_len = len(frame_number)
    new_frame_number = '0' * (digit_count - frame_number_len) + frame_number
    old_path = os.path.join(FRAMES_DIR, frame_name)
    new_path = os.path.join(FRAMES_DIR, 'out' + new_frame_number + '.png')
    os.rename(old_path, new_path)

