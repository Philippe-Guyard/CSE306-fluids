#bin/bash
rm frames/*
cmake . 
make 
./fluid_sim 
python3 rename_frames.py 3
ffmpeg -framerate 30 -pattern_type glob -i 'frames/*.png' -c:v libx264 -pix_fmt yuv420p out.mp4
